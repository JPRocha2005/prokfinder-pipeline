#!/usr/bin/env bash
#
# filter-low-quality-contigs-from-genomes.sh - Prokfinder module 3: contig filtering
#
# For each genome (seqkit):
#   - contigs shorter than <min_length> bp are removed
#   - the remaining contigs are sorted by length (longest first) and renamed Contig_1, Contig_2, ...
# The raw genomes are not changed; the clean genomes are written to another folder and are never moved from there
# (the next steps use their paths).
#
# Re-execution: a clean genome is reused if it is newer than its raw genome and was made with the same minimum
# length. Only new or changed genomes are processed.
#
# Usage: filter-low-quality-contigs-from-genomes.sh -l <accession_list> -i <raw_genomes_dir>
#                                                   -o <clean_genomes_dir> -t <table> [-m <min_length>] [-p <threads>] [-z]
#
# Output:
#   <clean_genomes_dir>/<ACCESSION>.fna[.gz]
#   <table>: ASSEMBLY_ACCESSION  CONTIGS_RAW  CONTIGS_CLEAN  MIN_CONTIG_LENGTH  CLEAN_GENOME_FILE  CONTIG_FILTER_STATUS
#            (status: OK | NO_CONTIGS_LEFT | FAILED_SEQKIT | NO_RAW_GENOME)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -i <raw_genomes_dir> -o <clean_genomes_dir> -t <table> [options]

  -l <file>   Accessions to process (one per line)
  -i <dir>    Folder with the raw genomes (<ACCESSION>.fna or .fna.gz)
  -o <dir>    Folder for the clean genomes
  -t <file>   Result table (also used as cache between runs)
  -m <bp>     Minimum contig length (default: 1000)
  -p <n>      Genomes processed in parallel (default: 4)
  -z          Write the clean genomes compressed (.fna.gz)
EOF
    exit "${1:-1}"
}

ACC_LIST=""
RAW_DIR=""
CLEAN_DIR=""
TABLE=""
MIN_LENGTH=1000
THREADS=4
GZIP=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -i) RAW_DIR="${2%/}"; shift 2 ;;
        -o) CLEAN_DIR="${2%/}"; shift 2 ;;
        -t) TABLE="$2"; shift 2 ;;
        -m) MIN_LENGTH="$2"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        -z) GZIP=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${ACC_LIST}" || -z "${RAW_DIR}" || -z "${CLEAN_DIR}" || -z "${TABLE}" ]]; then
    log "error" "Options -l, -i, -o and -t are required."
    usage
fi
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"
[[ "${MIN_LENGTH}" =~ ^[0-9]+$ ]] || die "The minimum contig length must be an integer (got '${MIN_LENGTH}')."
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "The number of threads must be a positive integer (got '${THREADS}')."

mkdir -p "${CLEAN_DIR}" "$(dirname "${TABLE}")"
setup_module_log "${CLEAN_DIR}/filter-low-quality-contigs-from-genomes.log"
require_programs seqkit

CLEAN_DIR="$(absolute_path "${CLEAN_DIR}")"
OUT_EXT="fna"
if [[ "${GZIP}" -eq 1 ]]; then
    OUT_EXT="fna.gz"
fi

HEADER=$'ASSEMBLY_ACCESSION\tCONTIGS_RAW\tCONTIGS_CLEAN\tMIN_CONTIG_LENGTH\tCLEAN_GENOME_FILE\tCONTIG_FILTER_STATUS'

make_temp_dir WORK_DIR
PENDING="${WORK_DIR}/pending.tsv"       # <accession> <TAB> <raw genome file>
NEW_ROWS="${WORK_DIR}/new_rows.tsv"
: > "${PENDING}"
: > "${NEW_ROWS}"

# Prints the number of sequences of a FASTA file (plain or .gz)
count_contigs() {
    local N
    if [[ "$1" == *.gz ]]; then
        N="$(gzip -dc -- "$1" | grep -c '^>' || true)"
    else
        N="$(grep -c '^>' -- "$1" || true)"
    fi
    echo "${N:-0}"
}

# Filters one genome. Runs in parallel (through xargs): a genome that fails is only recorded in the result file
# and the other genomes go on. The clean genome is written under a temporary name and renamed at the end, so an
# interrupted run never leaves a half-written genome that looks finished.
# Usage: filter_one_genome <clean_dir> <out_ext> <min_len> <result_file> <accession> <raw_file>
filter_one_genome() {
    local CLEAN_DIR="$1" OUT_EXT="$2" MIN_LEN="$3" RESULT_FILE="$4" ACC="$5" RAW="$6"
    local OUT="${CLEAN_DIR}/${ACC}.${OUT_EXT}"
    local TMP="${CLEAN_DIR}/${ACC}.tmp.${OUT_EXT}"
    local ERR BEFORE AFTER

    if ! ERR="$(set -o pipefail
                { seqkit seq --quiet --min-len "${MIN_LEN}" "${RAW}" \
                    | seqkit sort --quiet --by-length --reverse \
                    | seqkit replace --quiet --pattern '.*' --replacement 'Contig_{nr}' --out-file "${TMP}"; } 2>&1)"; then
        rm -f "${TMP}"
        printf '%s\tNULL\tNULL\t%s\tNULL\tFAILED_SEQKIT\n' "${ACC}" "${MIN_LEN}" >> "${RESULT_FILE}"
        log "warning" "seqkit failed on ${RAW}; the genome was skipped. $(printf '%s' "${ERR}" | tr '\n' ' ')"
        return 0
    fi

    BEFORE="$(count_contigs "${RAW}")"
    AFTER="$(count_contigs "${TMP}")"
    if [[ "${AFTER}" -eq 0 ]]; then
        rm -f "${TMP}" "${OUT}"
        printf '%s\t%s\t0\t%s\tNULL\tNO_CONTIGS_LEFT\n' "${ACC}" "${BEFORE}" "${MIN_LEN}" >> "${RESULT_FILE}"
        return 0
    fi

    mv -f "${TMP}" "${OUT}"
    # (one short line per genome: safe to append from several processes at once)
    printf '%s\t%s\t%s\t%s\t%s\tOK\n' "${ACC}" "${BEFORE}" "${AFTER}" "${MIN_LEN}" "${OUT}" >> "${RESULT_FILE}"
}

# 1. Genomes that can be reused from a previous run and genomes to process
declare -A CACHED_ROW=()
if [[ -s "${TABLE}" ]]; then
    while IFS=$'\t' read -r ACC REST; do
        CACHED_ROW["${ACC}"]="${REST}"
    done < <(tail -n +2 "${TABLE}")
fi

N_TOTAL=0
N_REUSED=0
while IFS= read -r ACC; do
    N_TOTAL=$((N_TOTAL + 1))
    RAW="$(genome_file_of "${RAW_DIR}" "${ACC}")"
    if [[ -z "${RAW}" ]]; then
        printf '%s\tNULL\tNULL\t%s\tNULL\tNO_RAW_GENOME\n' "${ACC}" "${MIN_LENGTH}" >> "${NEW_ROWS}"
        continue
    fi

    # Cached row: CONTIGS_RAW CONTIGS_CLEAN MIN_CONTIG_LENGTH CLEAN_GENOME_FILE STATUS
    ROW="${CACHED_ROW[${ACC}]:-}"
    if [[ -n "${ROW}" ]]; then
        IFS=$'\t' read -r _ _ CACHED_MIN CACHED_FILE CACHED_STATUS <<< "${ROW}"
        if [[ "${CACHED_MIN}" == "${MIN_LENGTH}" ]] && {
               { [[ "${CACHED_STATUS}" == "OK" && -f "${CACHED_FILE}" && "${CACHED_FILE}" -nt "${RAW}" ]]; } ||
               [[ "${CACHED_STATUS}" == "NO_CONTIGS_LEFT" ]]; }; then
            printf '%s\t%s\n' "${ACC}" "${ROW}" >> "${NEW_ROWS}"
            N_REUSED=$((N_REUSED + 1))
            continue
        fi
    fi
    printf '%s\t%s\n' "${ACC}" "${RAW}" >> "${PENDING}"
done < <(tr -d '\r' < "${ACC_LIST}" | awk 'NF { print $1 }' | sort -u)

N_PENDING="$(count_lines "${PENDING}")"
log "message" "Contig filtering (< ${MIN_LENGTH} bp removed): ${N_TOTAL} genome(s); ${N_REUSED} reused from a previous run; ${N_PENDING} to process (${THREADS} in parallel)."

# 2. Process the pending genomes in parallel
if [[ "${N_PENDING}" -gt 0 ]]; then
    log "message" "seqkit: $(seqkit version 2>&1 | head -n 1)"
    export -f log filter_one_genome count_contigs
    tr '\t' '\n' < "${PENDING}" \
        | xargs -d '\n' -n 2 -P "${THREADS}" bash -c 'filter_one_genome "$@"' _ \
            "${CLEAN_DIR}" "${OUT_EXT}" "${MIN_LENGTH}" "${NEW_ROWS}" \
        || die "Contig filtering stopped unexpectedly."
fi

# 3. Result table (rows of genomes that are not in this list are kept as cache for future runs)
{
    printf '%s\n' "${HEADER}"
    {
        cat "${NEW_ROWS}"
        if [[ -s "${TABLE}" ]]; then
            tail -n +2 "${TABLE}"
        fi
    } | awk -F'\t' '!seen[$1]++' | LC_ALL=C sort -t $'\t' -k1,1
} > "${TABLE}.tmp"
mv "${TABLE}.tmp" "${TABLE}"

N_OK="$(awk -F'\t' 'NR == FNR { want[$1] = 1; next } FNR > 1 && ($1 in want) && $6 == "OK"' \
        <(tr -d '\r' < "${ACC_LIST}") "${TABLE}" | wc -l)"
log "message" "Contig filtering complete: ${N_OK} of ${N_TOTAL} genome(s) have a clean genome in ${CLEAN_DIR}. Table: ${TABLE}"
