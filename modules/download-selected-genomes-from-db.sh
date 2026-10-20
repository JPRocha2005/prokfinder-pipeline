#!/usr/bin/env bash
#
# download-selected-genomes-from-db.sh - Prokfinder module 2: download of the selected genomes
#
# The NCBI reference download leaves one dehydrated manifest (fetch.txt) per domain, with the location of every
# genome file. A new manifest is built with the entries of the selected assemblies only, and
# "datasets rehydrate" downloads just those files. Only the genome sequence (*_genomic.fna) is kept: annotation
# files such as cds_from_genomic.fna or rna_from_genomic.fna are never downloaded.
# Each genome is saved as <ACCESSION>.fna (or .fna.gz with -z). Genomes already in the output folder are not
# downloaded again, so an interrupted run can be resumed.
#
# Usage: download-selected-genomes-from-db.sh -l <accession_list> -r <ncbi_dir> -o <raw_genomes_dir>
#                                             -t <status_table> [-d <domain_table>] [-j <threads>] [-z]
#
# Output:
#   <raw_genomes_dir>/<ACCESSION>.fna[.gz]
#   <status_table>: ASSEMBLY_ACCESSION  RAW_GENOME_FILE  DOWNLOAD_STATUS (DOWNLOADED | NOT_IN_FETCH | DOWNLOAD_FAILED)
#                   GENOME_DOMAIN_NCBI (BACTERIA | ARCHAEA: the fetch.txt that lists the genome)
#   <domain_table>: ASSEMBLY_ACCESSION  DOMAIN_NCBI (after two "#" lines that explain it)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -r <ncbi_dir> -o <raw_genomes_dir> -t <status_table> [-j <threads>] [-z]

  -l <file>   Accessions to download (one per line)
  -r <dir>    ~/prokfinder-refdb/NCBI (holds the fetch.txt manifests)
  -o <dir>    Folder where the genome files are saved
  -t <file>   Status table written at the end (one row per accession)
  -d <file>   Table of the NCBI domain (BACTERIA or ARCHAEA) of each genome, from the fetch.txt that lists it
  -j <n>      Concurrent downloads, 1-30 (default: 10)
  -z          Keep the genome files compressed (.fna.gz)

EOF
    exit "${1:-1}"
}

ACC_LIST=""
NCBI_DIR=""
OUTPUT_DIR=""
STATUS_TABLE=""
DOMAIN_TABLE=""
THREADS=10
GZIP=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -r) NCBI_DIR="${2%/}"; shift 2 ;;
        -o) OUTPUT_DIR="${2%/}"; shift 2 ;;
        -t) STATUS_TABLE="$2"; shift 2 ;;
        -d) DOMAIN_TABLE="$2"; shift 2 ;;
        -j) THREADS="$2"; shift 2 ;;
        -z) GZIP=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${ACC_LIST}" || -z "${NCBI_DIR}" || -z "${OUTPUT_DIR}" || -z "${STATUS_TABLE}" ]]; then
    log "error" "Options -l, -r, -o and -t are required."
    usage
fi
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"
if [[ ! "${THREADS}" =~ ^[0-9]+$ ]] || (( THREADS < 1 || THREADS > 30 )); then
    die "Invalid number of concurrent downloads '${THREADS}'. It must be an integer from 1 to 30."
fi

mkdir -p "${OUTPUT_DIR}"
setup_module_log "${OUTPUT_DIR}/download-selected-genomes-from-db.log"
require_programs datasets

REHYDRATE_MAX_ATTEMPTS=3    # attempts of "datasets rehydrate"
REHYDRATE_RETRY_WAIT=5      # seconds between attempts
REHYDRATION_DIR="${OUTPUT_DIR}/.rehydration"

make_temp_dir WORK_DIR
ALL_ACCESSIONS="${WORK_DIR}/all_accessions.txt"
PENDING_ACCESSIONS="${WORK_DIR}/pending_accessions.txt"
NOT_IN_FETCH="${WORK_DIR}/not_in_fetch.txt"
FAILED_DOWNLOAD="${WORK_DIR}/failed_download.txt"
: > "${NOT_IN_FETCH}"
: > "${FAILED_DOWNLOAD}"

# Runs "datasets rehydrate" on a directory that holds ncbi_dataset/fetch.txt. The command is repeated if it fails
# (e.g., transient network errors) up to REHYDRATE_MAX_ATTEMPTS times; files already downloaded are kept.
rehydrate() {
    local ARGS=(--no-progressbar --max-workers "${THREADS}" --directory "${REHYDRATION_DIR}")
    local ATTEMPT=1

    if [[ "${GZIP}" -eq 1 ]]; then
        ARGS+=(--gzip)
    fi

    until datasets rehydrate "${ARGS[@]}"; do
        if (( ATTEMPT >= REHYDRATE_MAX_ATTEMPTS )); then
            log "warning" "datasets rehydrate failed ${ATTEMPT} time(s). The genomes not downloaded are retried in the next run."
            return 1
        fi
        log "message" "datasets rehydrate failed or was interrupted (attempt ${ATTEMPT}/${REHYDRATE_MAX_ATTEMPTS}); retrying in ${REHYDRATE_RETRY_WAIT}s..."
        ATTEMPT=$((ATTEMPT + 1))
        sleep "${REHYDRATE_RETRY_WAIT}"
    done
}

# 1. Accessions requested and accessions still missing in the output folder
tr -d '\r' < "${ACC_LIST}" | awk 'NF { print $1 }' | sort -u > "${ALL_ACCESSIONS}"
while IFS= read -r ACC; do
    if [[ -z "$(genome_file_of "${OUTPUT_DIR}" "${ACC}")" ]]; then
        echo "${ACC}"
    fi
done < "${ALL_ACCESSIONS}" > "${PENDING_ACCESSIONS}"

TOTAL_COUNT="$(count_lines "${ALL_ACCESSIONS}")"
PENDING_COUNT="$(count_lines "${PENDING_ACCESSIONS}")"
log "message" "Genome download: ${TOTAL_COUNT} selected; $((TOTAL_COUNT - PENDING_COUNT)) already in ${OUTPUT_DIR}; ${PENDING_COUNT} to download."

# 2. One pass over the fetch.txt manifests (one per NCBI domain) finds, for every selected accession, the domain
#    whose manifest lists it (BACTERIA or ARCHAEA, by the NCBI taxonomy) and, for the pending accessions, the line
#    of the genome file. A fetch.txt line has 3 tab-separated columns: an opaque URL, a number and the file path
#    (data/<accession>/<file>). The accession is read from the path and compared exactly (GCA_000000001.1 must not
#    match GCA_000000001.15). Only the genome sequence is kept. The manifests written by NCBI end their lines with
#    CRLF: the \r is removed before the path is read.
FETCH_FILES=()
FETCH_DOMAINS=""
for DOMAIN in "${PF_DOMAINS[@]}"; do
    FETCH_FILE="$(pf_ncbi_fetch_file "${NCBI_DIR}" "${DOMAIN}")"
    if [[ -s "${FETCH_FILE}" ]]; then
        FETCH_FILES+=("${FETCH_FILE}")
        FETCH_DOMAINS="${FETCH_DOMAINS:+${FETCH_DOMAINS},}${DOMAIN}"
    fi
done
if [[ "${#FETCH_FILES[@]}" -eq 0 ]]; then
    die "No fetch.txt manifest found under ${NCBI_DIR}. Was the NCBI reference database downloaded?"
fi

rm -rf "${REHYDRATION_DIR}"
mkdir -p "${REHYDRATION_DIR}/ncbi_dataset/data"
MANIFEST="${REHYDRATION_DIR}/ncbi_dataset/fetch.txt"
DOMAIN_OF="${WORK_DIR}/domain_of_accession.tsv"

log "message" "Looking for the selected genomes in the fetch.txt manifests (${FETCH_DOMAINS//,/ and })..."
awk -F'\t' -v all_file="${ALL_ACCESSIONS}" -v pending_file="${PENDING_ACCESSIONS}" \
    -v missing_file="${NOT_IN_FETCH}" -v domain_file="${DOMAIN_OF}" -v domains="${FETCH_DOMAINS}" '
    BEGIN {
        while ((getline a < all_file) > 0) if (a != "") selected[a] = 1
        close(all_file)
        while ((getline a < pending_file) > 0) if (a != "") pending[a] = 1
        close(pending_file)
        split(domains, domain_of_file, ",")
    }
    FNR == 1 { file_index++ }
    {
        sub(/\r$/, "")
        path = (NF >= 3 ? $3 : $0)
        if (path !~ /_genomic[.]fna$/ || path ~ /(cds|rna)_from_genomic[.]fna$/) next
        n = split(path, parts, "/")
        acc = (n >= 2 ? parts[n - 1] : "")
        if (!(acc in selected) || (acc in found)) next
        found[acc] = 1
        print acc "\t" domain_of_file[file_index] > domain_file
        if (acc in pending) print
    }
    END {
        for (a in selected) if (!(a in found)) print a > missing_file
        close(domain_file)
    }
' "${FETCH_FILES[@]}" > "${MANIFEST}"
touch "${DOMAIN_OF}"

if [[ -s "${NOT_IN_FETCH}" ]]; then
    log "warning" "$(count_lines "${NOT_IN_FETCH}") selected accession(s) have no genome file in the fetch.txt manifests (${FETCH_DOMAINS//,/ and }); the ones not downloaded yet cannot be downloaded."
fi

if [[ "${PENDING_COUNT}" -gt 0 ]]; then
    log "message" "Rehydration manifest prepared with $(count_lines "${MANIFEST}") genome file(s)."

    # 3. Download, then move each genome to <OUTPUT_DIR>/<ACCESSION>.fna[.gz]
    if [[ -s "${MANIFEST}" ]]; then
        log "message" "Rehydrating genomes ($(datasets --version 2>&1 | head -n 1); ${THREADS} concurrent downloads)..."
        rehydrate || true

        MOVED=0
        while IFS= read -r FILE; do
            ACC="$(basename "$(dirname "${FILE}")")"
            if [[ "${FILE}" == *.gz ]]; then
                mv -- "${FILE}" "${OUTPUT_DIR}/${ACC}.fna.gz"
            else
                mv -- "${FILE}" "${OUTPUT_DIR}/${ACC}.fna"
            fi
            MOVED=$((MOVED + 1))
        done < <(find "${REHYDRATION_DIR}/ncbi_dataset/data" -mindepth 2 -maxdepth 2 -type f \
                    \( -name '*_genomic.fna' -o -name '*_genomic.fna.gz' \) \
                    ! -name '*_from_genomic.fna*')
        log "message" "${MOVED} genome file(s) downloaded."
    fi
fi
rm -rf "${REHYDRATION_DIR}"

# 4. Status and NCBI domain of every requested accession
declare -A DOMAIN_OF_ACC=()
while IFS=$'\t' read -r ACC DOMAIN; do
    DOMAIN_OF_ACC["${ACC}"]="${DOMAIN}"
done < "${DOMAIN_OF}"
{
    printf 'ASSEMBLY_ACCESSION\tRAW_GENOME_FILE\tDOWNLOAD_STATUS\tGENOME_DOMAIN_NCBI\n'
    while IFS= read -r ACC; do
        FILE="$(genome_file_of "${OUTPUT_DIR}" "${ACC}")"
        DOMAIN="${DOMAIN_OF_ACC[${ACC}]:-NULL}"
        if [[ -n "${FILE}" ]]; then
            printf '%s\t%s\tDOWNLOADED\t%s\n' "${ACC}" "$(absolute_path "${FILE}")" "${DOMAIN}"
        elif grep -qxF "${ACC}" "${NOT_IN_FETCH}"; then
            printf '%s\tNULL\tNOT_IN_FETCH\t%s\n' "${ACC}" "${DOMAIN}"
        else
            printf '%s\tNULL\tDOWNLOAD_FAILED\t%s\n' "${ACC}" "${DOMAIN}"
        fi
    done < "${ALL_ACCESSIONS}"
} > "${STATUS_TABLE}"

# Domain of each genome by the NCBI taxonomy (used to choose the domain-specific tools of the next steps)
if [[ -n "${DOMAIN_TABLE}" ]]; then
    mkdir -p "$(dirname "${DOMAIN_TABLE}")"
    {
        echo "# Domain of each genome according to the NCBI taxonomy (the NCBI fetch.txt manifest - BACTERIA or ARCHAEA - that lists it)."
        echo "# The GTDB taxonomy may place a genome in another domain."
        printf 'ASSEMBLY_ACCESSION\tDOMAIN_NCBI\n'
        sort -k1,1 "${DOMAIN_OF}"
    } > "${DOMAIN_TABLE}"
    log "message" "NCBI domain of the genomes: $(awk -F'\t' '$2 == "BACTERIA"' "${DOMAIN_OF}" | wc -l) BACTERIA, $(awk -F'\t' '$2 == "ARCHAEA"' "${DOMAIN_OF}" | wc -l) ARCHAEA -> ${DOMAIN_TABLE}"
fi

N_OK="$(awk -F'\t' 'NR > 1 && $3 == "DOWNLOADED"' "${STATUS_TABLE}" | wc -l)"
N_FAILED="$(awk -F'\t' 'NR > 1 && $3 == "DOWNLOAD_FAILED"' "${STATUS_TABLE}" | wc -l)"
log "message" "Genome files available: ${N_OK} of ${TOTAL_COUNT}. Status table: ${STATUS_TABLE}"
if [[ "${N_FAILED}" -gt 0 ]]; then
    log "warning" "${N_FAILED} genome(s) could not be downloaded in this run; run Prokfinder again to retry them."
fi
