#!/usr/bin/env bash
#
# download-reference-db-gtdb.sh - Prokfinder module 0.3: GTDB reference database
#
# Downloads the taxonomy and the metadata of the genomes classified by GTDB (latest release by default):
#     ar53_taxonomy.tsv.gz   ar53_metadata.tsv.gz   bac120_taxonomy.tsv.gz   bac120_metadata.tsv.gz
# and formats them:
#   - taxonomy: the lineage "d__...;p__...;...;s__..." is split into one column per rank (prefixes removed);
#     ranks with no name (e.g. "s__") -> NULL
#   - metadata: header without spaces, empty/invalid cells -> NULL; the GTDB column "accession" (with the RS_ / GB_
#     prefix) is renamed gtdb_accession, so it is not confused with ACCESSION
#   - both: a first column ACCESSION with the NCBI assembly accession (GTDB adds the prefixes RS_ / GB_)
#
# Tables written by an older Prokfinder (older format in the checkpoint) are rebuilt from the .tsv.gz files already
# downloaded.
#
# Usage: download-reference-db-gtdb.sh -o <gtdb_dir> [--url <release_url>] [--force]
#
# Output inside <gtdb_dir>:
#   ARCHAEA-GTDB-TAXONOMY/ar53_taxonomy.tsv.gz + ar53_taxonomy_formatted.tsv
#   BACTERIA-GTDB-TAXONOMY/bac120_taxonomy.tsv.gz + bac120_taxonomy_formatted.tsv
#   ARCHAEA-GTDB-METADATA/ar53_metadata.tsv.gz + ar53_metadata_formatted.tsv
#   BACTERIA-GTDB-METADATA/bac120_metadata.tsv.gz + bac120_metadata_formatted.tsv
#   VERSION.txt              (GTDB release)
#   COMPLETED-DOWNLOAD.txt   (checkpoint: written only at the end)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -o <gtdb_dir> [--url <release_url>] [--force]

  -o <gtdb_dir>     Directory where the GTDB reference data will be saved (required)
  --url <url>       GTDB release folder (default: ${GTDB_DEFAULT_URL})
  --force           Remove the previous download of this module and start again
  -h                Show this help message
EOF
    exit "${1:-1}"
}

GTDB_DEFAULT_URL="https://data.gtdb.aau.ecogenomic.org/releases/latest"
GTDB_URL="${GTDB_DEFAULT_URL}"
OUTPUT_DIR=""
FORCE=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -o) OUTPUT_DIR="${2%/}"; shift 2 ;;
        --url) GTDB_URL="${2%/}"; shift 2 ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${OUTPUT_DIR}" ]]; then
    log "error" "Output directory (-o) is required."
    usage
fi

mkdir -p "${OUTPUT_DIR}"
setup_module_log "${OUTPUT_DIR}/download-reference-db-gtdb.log"
require_programs curl gzip md5sum

COMPLETED_MARKER="$(pf_completed_marker "${OUTPUT_DIR}")"
RANKS="DOMAIN PHYLUM CLASS ORDER FAMILY GENUS SPECIES"

############ FUNCTIONS ################

# Downloads <url> to <file> (resumes partial downloads and retries on network errors)
download_file() {
    local URL="$1"
    local FILE="$2"
    log "message" "Downloading ${URL}"
    curl -fSL --retry 5 --retry-delay 10 -C - -o "${FILE}.part" "${URL}" \
        || die "Could not download ${URL}"
    mv "${FILE}.part" "${FILE}"
}

# Checks <file> against the MD5SUM.txt of the release (skipped if the release has no checksum for it)
check_md5() {
    local FILE="$1"
    local SUMS="$2"
    local EXPECTED
    # MD5SUM.txt lists the versioned names (e.g. ar53_taxonomy_r232.tsv.gz), while "latest" serves ar53_taxonomy.tsv.gz
    EXPECTED="$(awk -v f="$(basename "${FILE}")" '
        { n = $2; sub(/^.*\//, "", n); m = n; sub(/_r[0-9]+[.]/, ".", m); if (n == f || m == f) { print $1; exit } }
    ' "${SUMS}" 2>/dev/null || true)"
    if [[ -z "${EXPECTED}" ]]; then
        log "warning" "No checksum found for $(basename "${FILE}") in MD5SUM.txt; integrity not verified."
        return 0
    fi
    if [[ "$(md5sum "${FILE}" | cut -d' ' -f1)" != "${EXPECTED}" ]]; then
        rm -f "${FILE}"
        die "Checksum mismatch for ${FILE} (the file was removed; run the module again)."
    fi
}

# Taxonomy table: <gtdb_accession> <TAB> d__X;p__Y;...;s__Z  ->  formatted TSV
format_taxonomy() {
    gzip -dc "$1" | awk -F'\t' -v OFS='\t' -v ranks="${RANKS}" '
    BEGIN {
        nr = split(ranks, rank, " ")
        out = "ACCESSION" OFS "GTDB_ACCESSION"
        for (i = 1; i <= nr; i++) out = out OFS "GTDB_" rank[i]
        print out, "GTDB_TAXONOMY"
    }
    {
        gsub(/\r/, "")
        acc = $1
        sub(/^(RS|GB)_/, "", acc)
        out = acc OFS $1
        n = split($2, lin, ";")
        for (i = 1; i <= nr; i++) {
            v = (i <= n ? lin[i] : "")
            sub(/^[a-z]__/, "", v)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
            out = out OFS (v == "" ? "NULL" : v)
        }
        print out, $2
    }'
}

# Metadata table: header kept (spaces -> "_"; "accession" -> "gtdb_accession"), ACCESSION added as first column,
# invalid cells -> NULL
format_metadata() {
    gzip -dc "$1" | format_tsv 0 | awk -F'\t' -v OFS='\t' '
        NR == 1 { if (tolower($1) == "accession") $1 = "gtdb_accession"; print "ACCESSION", $0; next }
        { acc = $1; sub(/^(RS|GB)_/, "", acc); print acc, $0 }'
}

############ DOWNLOAD ################

if [[ "${FORCE}" -eq 1 ]]; then
    log "message" "Option --force: removing the previous GTDB reference data in ${OUTPUT_DIR}"
    rm -f "${COMPLETED_MARKER}" "$(pf_gtdb_version_file "${OUTPUT_DIR}")"
    for DOMAIN in "${PF_DOMAINS[@]}"; do
        rm -rf "${OUTPUT_DIR}/${DOMAIN}-GTDB-TAXONOMY" "${OUTPUT_DIR}/${DOMAIN}-GTDB-METADATA"
    done
elif [[ -f "${COMPLETED_MARKER}" ]] && ! pf_refdb_format_ok "${OUTPUT_DIR}" "${PF_REFDB_FORMAT_GTDB}"; then
    # Tables of an older Prokfinder: they are formatted again from the files already downloaded
    log "message" "The GTDB tables were written by an older version of Prokfinder: formatting them again from the files already downloaded."
    for DOMAIN in "${PF_DOMAINS[@]}"; do
        rm -f "$(pf_gtdb_taxonomy_table "${OUTPUT_DIR}" "${DOMAIN}")" "$(pf_gtdb_metadata_table "${OUTPUT_DIR}" "${DOMAIN}")"
    done
fi

log "message" "Downloading GTDB reference data from ${GTDB_URL} into: ${OUTPUT_DIR}"

VERSION_FILE="$(pf_gtdb_version_file "${OUTPUT_DIR}")"
SUMS_FILE="${OUTPUT_DIR}/MD5SUM.txt"

# The release is read again only when a file has to be downloaded (otherwise the release of the files on disk is kept)
MISSING_RAW=0
for DOMAIN in "${PF_DOMAINS[@]}"; do
    PREFIX="$(pf_gtdb_prefix "${DOMAIN}")"
    [[ -s "${OUTPUT_DIR}/${DOMAIN}-GTDB-TAXONOMY/${PREFIX}_taxonomy.tsv.gz" ]] || MISSING_RAW=1
    [[ -s "${OUTPUT_DIR}/${DOMAIN}-GTDB-METADATA/${PREFIX}_metadata.tsv.gz" ]] || MISSING_RAW=1
done
if [[ "${MISSING_RAW}" -eq 1 || ! -s "${VERSION_FILE}" ]]; then
    curl -fsSL --retry 3 -o "${VERSION_FILE}" "${GTDB_URL}/VERSION.txt" || echo "unknown" > "${VERSION_FILE}"
    curl -fsSL --retry 3 -o "${SUMS_FILE}" "${GTDB_URL}/MD5SUM.txt" || : > "${SUMS_FILE}"
else
    : > "${SUMS_FILE}"
fi
log "message" "GTDB release: $(head -n 1 "${VERSION_FILE}")"

for DOMAIN in "${PF_DOMAINS[@]}"; do
    PREFIX="$(pf_gtdb_prefix "${DOMAIN}")"

    for KIND in taxonomy metadata; do
        if [[ "${KIND}" == "taxonomy" ]]; then
            DIR="${OUTPUT_DIR}/${DOMAIN}-GTDB-TAXONOMY"
            FORMATTED="$(pf_gtdb_taxonomy_table "${OUTPUT_DIR}" "${DOMAIN}")"
        else
            DIR="${OUTPUT_DIR}/${DOMAIN}-GTDB-METADATA"
            FORMATTED="$(pf_gtdb_metadata_table "${OUTPUT_DIR}" "${DOMAIN}")"
        fi
        RAW="${DIR}/${PREFIX}_${KIND}.tsv.gz"

        if [[ -s "${FORMATTED}" ]]; then
            log "message" "[${DOMAIN}] ${KIND}: already formatted. Skipping..."
            continue
        fi

        mkdir -p "${DIR}"
        if [[ ! -s "${RAW}" ]]; then
            download_file "${GTDB_URL}/${PREFIX}_${KIND}.tsv.gz" "${RAW}"
            check_md5 "${RAW}" "${SUMS_FILE}"
        fi

        log "message" "[${DOMAIN}] Formatting ${KIND} table..."
        if [[ "${KIND}" == "taxonomy" ]]; then
            format_taxonomy "${RAW}" > "${FORMATTED}.tmp"
        else
            format_metadata "${RAW}" > "${FORMATTED}.tmp"
        fi
        mv "${FORMATTED}.tmp" "${FORMATTED}"
        log "message" "[${DOMAIN}] ${KIND}: $(count_records "${FORMATTED}") genomes -> ${FORMATTED}"
    done
done

rm -f "${SUMS_FILE}"

{
    echo "GTDB reference data downloaded on $(date +"%Y-%m-%d %H:%M:%S")"
    echo "url=${GTDB_URL}"
    echo "release=$(head -n 1 "${VERSION_FILE}")"
    echo "format=${PF_REFDB_FORMAT_GTDB}"
} > "${COMPLETED_MARKER}"

log "message" "GTDB reference data download completed successfully!"
