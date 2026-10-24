#!/usr/bin/env bash
#
# download-reference-db-gtdb.sh - Prokfinder module 0.3: GTDB reference database
#
# Downloads the taxonomy and the metadata of the genomes classified by GTDB (latest release by default):
#     ar53_taxonomy.tsv.gz   ar53_metadata.tsv.gz   bac120_taxonomy.tsv.gz   bac120_metadata.tsv.gz
# and formats them:
#   - taxonomy: the lineage "d__...;p__...;...;s__..." is split into one column per rank (prefixes removed);
#     ranks with no name (e.g. "s__") -> NULL. Four columns of the metadata are added, on the status of the genome
#     in GTDB: GTDB_SPECIES_REPRESENTATIVE (yes/no: the genome is the representative of its GTDB species),
#     GTDB_REPRESENTATIVE_GENOME (accession of that representative), GTDB_TYPE_MATERIAL (e.g. "type strain of
#     species"; NULL when it is not type material) and GTDB_TYPE_SPECIES_OF_GENUS (yes/no)
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

# Status of each genome in GTDB, from the raw metadata table: <gtdb_accession> <TAB> representative (yes/no) <TAB>
# accession of the species representative <TAB> type material <TAB> type species of genus (yes/no)
# Empty or invalid values (e.g., "not type material", "not used as type") become NULL with the rule of every
# reference table (PF_AWK_NULL); the converted values are appended to <null_rows>
# Usage: gtdb_status <metadata.tsv.gz> <null_rows> <table name>
gtdb_status() {
    gzip -dc "$1" | awk -F'\t' -v OFS='\t' -v null_rows="$2" -v table="$3" "${PF_AWK_NULL}"'
    function yes_no(x) { x = tolower(x); return (x == "t" || x == "true") ? "yes" : ((x == "f" || x == "false") ? "no" : "NULL") }
    function value(x, column,   rule) {
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", x)
        if ((rule = pf_null_rule(x)) == "") return x
        pf_null_record(column, x, rule)
        return "NULL"
    }
    END { pf_null_write(null_rows, table) }
    NR == 1 {
        gsub(/\r/, "")
        for (i = 1; i <= NF; i++) col[$i] = i
        if (!("accession" in col) || !("gtdb_representative" in col) || !("gtdb_genome_representative" in col) ||
            !("gtdb_type_designation_ncbi_taxa" in col) || !("gtdb_type_species_of_genus" in col)) {
            print "GTDB metadata: the columns accession, gtdb_representative, gtdb_genome_representative, gtdb_type_designation_ncbi_taxa and gtdb_type_species_of_genus are required" > "/dev/stderr"
            exit 2
        }
        next
    }
    {
        gsub(/\r/, "")
        rep = $(col["gtdb_genome_representative"]); sub(/^(RS|GB)_/, "", rep)
        print $(col["accession"]), yes_no($(col["gtdb_representative"])), value(rep, "GTDB_REPRESENTATIVE_GENOME"),
              value($(col["gtdb_type_designation_ncbi_taxa"]), "GTDB_TYPE_MATERIAL"), yes_no($(col["gtdb_type_species_of_genus"]))
    }'
}

# Taxonomy table: <gtdb_accession> <TAB> d__X;p__Y;...;s__Z  ->  formatted TSV, with the status of each genome in
# GTDB (from <status>, written by gtdb_status)
# (a rank without name, e.g. "s__", becomes NULL and is appended to <null_rows> as an EMPTY value)
# Usage: format_taxonomy <taxonomy.tsv.gz> <status> <null_rows> <table name>
format_taxonomy() {
    gzip -dc "$1" | awk -F'\t' -v OFS='\t' -v ranks="${RANKS}" -v status="$2" -v null_rows="$3" -v table="$4" "${PF_AWK_NULL}"'
    BEGIN {
        nr = split(ranks, rank, " ")
        while ((getline line < status) > 0) { split(line, f, "\t"); st[f[1]] = f[2] OFS f[3] OFS f[4] OFS f[5] }
        close(status)
        out = "ACCESSION" OFS "GTDB_ACCESSION"
        for (i = 1; i <= nr; i++) out = out OFS "GTDB_" rank[i]
        print out, "GTDB_TAXONOMY", "GTDB_SPECIES_REPRESENTATIVE", "GTDB_REPRESENTATIVE_GENOME", "GTDB_TYPE_MATERIAL",
              "GTDB_TYPE_SPECIES_OF_GENUS"
    }
    {
        gsub(/\r/, "")
        acc = $1
        sub(/^(RS|GB)_/, "", acc)
        out = acc OFS $1
        n = split($2, lin, ";")
        for (i = 1; i <= nr; i++) {
            v = (i <= n ? lin[i] : "")
            raw = v
            sub(/^[a-z]__/, "", v)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
            if ((rule = pf_null_rule(v)) != "") { pf_null_record("GTDB_" rank[i], raw, rule); v = "NULL" }
            out = out OFS v
        }
        print out, $2, (($1 in st) ? st[$1] : "NULL" OFS "NULL" OFS "NULL" OFS "NULL")
    }
    END { pf_null_write(null_rows, table) }'
}

# Metadata table: header kept (spaces -> "_"; "accession" -> "gtdb_accession"), ACCESSION added as first column,
# invalid cells -> NULL (appended to <null_rows>)
# Usage: format_metadata <metadata.tsv.gz> <null_rows> <table name>
format_metadata() {
    gzip -dc "$1" | format_tsv 0 "$2" "$3" | awk -F'\t' -v OFS='\t' '
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

    # (the metadata first: the taxonomy table takes the status of each genome from it)
    for KIND in metadata taxonomy; do
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
        # (values converted to NULL: one part per table, kept next to it and joined in NULL-VALUES-GTDB.tsv)
        NULL_ROWS="${FORMATTED%.tsv}.null-values.part.tsv"
        TABLE_NAME="$(basename "${DIR}")/$(basename "${FORMATTED}")"
        rm -f "${NULL_ROWS}"
        if [[ "${KIND}" == "taxonomy" ]]; then
            STATUS="${FORMATTED}.status.tmp"
            gtdb_status "${OUTPUT_DIR}/${DOMAIN}-GTDB-METADATA/${PREFIX}_metadata.tsv.gz" "${NULL_ROWS}" "${TABLE_NAME}" \
                > "${STATUS}" || die "Could not read the GTDB metadata of ${DOMAIN} (its columns may have changed)."
            format_taxonomy "${RAW}" "${STATUS}" "${NULL_ROWS}" "${TABLE_NAME}" > "${FORMATTED}.tmp"
            rm -f "${STATUS}"
        else
            format_metadata "${RAW}" "${NULL_ROWS}" "${TABLE_NAME}" > "${FORMATTED}.tmp"
        fi
        mv "${FORMATTED}.tmp" "${FORMATTED}"
        log "message" "[${DOMAIN}] ${KIND}: $(count_records "${FORMATTED}") genomes -> ${FORMATTED}"
    done
done

rm -f "${SUMS_FILE}"

pf_write_null_mapping "${OUTPUT_DIR}/NULL-VALUES-GTDB.tsv" "${OUTPUT_DIR}"/*-GTDB-*/*.null-values.part.tsv
log "message" "Values converted to NULL (empty or invalid): ${OUTPUT_DIR}/NULL-VALUES-GTDB.tsv"

{
    echo "GTDB reference data downloaded on $(date +"%Y-%m-%d %H:%M:%S")"
    echo "url=${GTDB_URL}"
    echo "release=$(head -n 1 "${VERSION_FILE}")"
    echo "format=${PF_REFDB_FORMAT_GTDB}"
} > "${COMPLETED_MARKER}"

log "message" "GTDB reference data download completed successfully!"
