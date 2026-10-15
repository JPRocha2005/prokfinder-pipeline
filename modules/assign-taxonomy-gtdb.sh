#!/usr/bin/env bash
#
# assign-taxonomy-gtdb.sh - Prokfinder module 5: GTDB taxonomy of the genomes
#
# 1. The taxonomy of each genome is looked up in the GTDB taxonomy tables (phylogenomic tree of the database).
#    GTDB lists each genome under only one of its accessions, so the lookup tries, in this order:
#    the same accession; the paired GenBank/RefSeq accession (GCA_x.1 <-> GCF_x.1); another version (GCA_x.*).
# 2. Genomes not found in GTDB:
#    - by default they are written to NON-CLASSIFIED-GENOMES-BY-GTDB.txt and left out of the next steps
#    - with --gtdbtk they are classified with GTDB-Tk (classify_wf: placement in the GTDB tree).
#      WARNING: GTDB-Tk needs a large amount of RAM (> 200 GB with the full reference package)
#
# Re-execution: GTDB-Tk results are kept in a cache table, so each genome is classified by GTDB-Tk only once.
#
# Usage: assign-taxonomy-gtdb.sh -l <accession_list> -g <gtdb_dir> -o <taxonomy_dir>
#                                [--gtdbtk -i <clean_genomes_dir> -p <threads>]
#
# Output inside <taxonomy_dir>:
#   genome-taxonomy.tsv          ASSEMBLY_ACCESSION GTDB_DOMAIN ... GTDB_SPECIES GTDB_TAXONOMY_SOURCE
#   NON-CLASSIFIED-GENOMES-BY-GTDB.txt   accessions without GTDB taxonomy
#   GTDBTK-RESULTS/              GTDB-Tk runs and gtdbtk-taxonomy-cache.tsv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -g <gtdb_dir> -o <taxonomy_dir> [--gtdbtk -i <clean_genomes_dir>] [-p <threads>]

  -l <file>    Accessions to classify (one per line)
  -g <dir>     ~/prokfinder-refdb/GTDB
  -o <dir>     Output folder
  --gtdbtk     Classify the genomes not found in the GTDB tables with GTDB-Tk (needs > 200 GB of RAM)
  -i <dir>     Folder with the clean genomes (required with --gtdbtk)
  -p <n>       CPUs used by GTDB-Tk (default: 4)

Environment variables: GTDBTK_DATA_PATH (required by GTDB-Tk), GTDBTK_EXTRA_ARGS (extra options for classify_wf,
e.g. "--scratch_dir <dir>" to reduce the pplacer memory; GTDB-Tk 2.7 needs no extra option)
EOF
    exit "${1:-1}"
}

ACC_LIST=""
GTDB_DIR=""
TAXONOMY_DIR=""
CLEAN_DIR=""
USE_GTDBTK=0
THREADS=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -g) GTDB_DIR="${2%/}"; shift 2 ;;
        -o) TAXONOMY_DIR="${2%/}"; shift 2 ;;
        -i) CLEAN_DIR="${2%/}"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        --gtdbtk) USE_GTDBTK=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${ACC_LIST}" || -z "${GTDB_DIR}" || -z "${TAXONOMY_DIR}" ]]; then
    log "error" "Options -l, -g and -o are required."
    usage
fi
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"
if [[ "${USE_GTDBTK}" -eq 1 && -z "${CLEAN_DIR}" ]]; then
    die "Option -i (clean genomes folder) is required with --gtdbtk."
fi

mkdir -p "${TAXONOMY_DIR}"
setup_module_log "${TAXONOMY_DIR}/assign-taxonomy-gtdb.log"

TAXONOMY_TABLE="${TAXONOMY_DIR}/genome-taxonomy.tsv"
NON_CLASSIFIED_LIST="${TAXONOMY_DIR}/NON-CLASSIFIED-GENOMES-BY-GTDB.txt"
GTDBTK_DIR="${TAXONOMY_DIR}/GTDBTK-RESULTS"
GTDBTK_CACHE="${GTDBTK_DIR}/gtdbtk-taxonomy-cache.tsv"
RANKS="DOMAIN PHYLUM CLASS ORDER FAMILY GENUS SPECIES"
GTDB_RELEASE="$(head -n 1 "$(pf_gtdb_version_file "${GTDB_DIR}")" 2>/dev/null || echo "unknown")"

GTDB_TABLES=()
for DOMAIN in "${PF_DOMAINS[@]}"; do
    TABLE="$(pf_gtdb_taxonomy_table "${GTDB_DIR}" "${DOMAIN}")"
    [[ -f "${TABLE}" ]] || die "GTDB taxonomy table not found: ${TABLE}. Was the GTDB reference database downloaded?"
    GTDB_TABLES+=("${TABLE}")
done

make_temp_dir WORK_DIR
tr -d '\r' < "${ACC_LIST}" | awk 'NF { print $1 }' | sort -u > "${WORK_DIR}/accessions.txt"
log "message" "Taxonomic classification of $(count_lines "${WORK_DIR}/accessions.txt") genome(s) with the GTDB taxonomy (release ${GTDB_RELEASE})..."

# 1. Lookup in the GTDB tables. Only the rows of the requested genomes are kept in memory.
#    Output: <accession> <TAB> 7 ranks <TAB> source
awk -F'\t' -v OFS='\t' -v acc_file="${WORK_DIR}/accessions.txt" -v release="${GTDB_RELEASE}" -v ranks="${RANKS}" '
    function num(a) { sub(/^GC[AF]_/, "", a); return a }
    function base(a) { a = num(a); sub(/[.][0-9]+$/, "", a); return a }
    BEGIN {
        nr = split(ranks, rank, " ")
        while ((getline a < acc_file) > 0) { want[++n] = a; by_exact[a] = 1; by_num[num(a)] = 1; by_base[base(a)] = 1 }
        close(acc_file)
    }
    FNR == 1 {
        for (i = 1; i <= NF; i++) for (k = 1; k <= nr; k++) if ($i == "GTDB_" rank[k]) col[k] = i
        next
    }
    ($1 in by_exact) || (num($1) in by_num) || (base($1) in by_base) {
        line = ""
        for (k = 1; k <= nr; k++) line = line OFS $(col[k])
        exact[$1] = line; paired[num($1)] = line; other[base($1)] = line
    }
    END {
        out = "ASSEMBLY_ACCESSION"
        for (k = 1; k <= nr; k++) out = out OFS "GTDB_" rank[k]
        print out, "GTDB_TAXONOMY_SOURCE"
        for (i = 1; i <= n; i++) {
            a = want[i]
            if (a in exact)              print a exact[a], "GTDB_" release
            else if (num(a) in paired)   print a paired[num(a)], "GTDB_" release "_PAIRED_ACCESSION"
            else if (base(a) in other)   print a other[base(a)], "GTDB_" release "_OTHER_VERSION"
            else {
                line = a
                for (k = 1; k <= nr; k++) line = line OFS "NULL"
                print line, "NULL"
            }
        }
    }
' "${GTDB_TABLES[@]}" > "${WORK_DIR}/taxonomy.tsv"

awk -F'\t' 'NR > 1 && $NF == "NULL" { print $1 }' "${WORK_DIR}/taxonomy.tsv" > "${WORK_DIR}/not_in_gtdb.txt"
N_NOT_IN_GTDB="$(count_lines "${WORK_DIR}/not_in_gtdb.txt")"
log "message" "Found in the GTDB tables: $(( $(count_records "${WORK_DIR}/taxonomy.tsv") - N_NOT_IN_GTDB )); not found: ${N_NOT_IN_GTDB}"

# 2. GTDB-Tk for the genomes not found in the tables (optional)
#    Parses the summary files of a GTDB-Tk run into "<accession> <TAB> 7 ranks <TAB> source"
parse_gtdbtk_summaries() {
    local RUN_DIR="$1"
    local SOURCE="$2"
    find "${RUN_DIR}" -maxdepth 2 -type f -name 'gtdbtk.*.summary.tsv' -print0 \
        | xargs -0 -r awk -F'\t' -v OFS='\t' -v source="${SOURCE}" '
            FNR == 1 { for (i = 1; i <= NF; i++) { if ($i == "user_genome") cg = i; if ($i == "classification") cc = i }; next }
            {
                if (!match($cg, /GC[AF]_[0-9]+[.][0-9]+/)) next
                acc = substr($cg, RSTART, RLENGTH)
                if (acc in done) next
                done[acc] = 1
                if ($cc ~ /^Unclassified/ || $cc == "") next
                n = split($cc, lin, ";")
                out = acc
                for (k = 1; k <= 7; k++) { v = (k <= n ? lin[k] : ""); sub(/^[a-z]__/, "", v); out = out OFS (v == "" ? "NULL" : v) }
                print out, source
            }'
}

if [[ "${USE_GTDBTK}" -eq 1 && "${N_NOT_IN_GTDB}" -gt 0 ]]; then
    require_programs gtdbtk
    [[ -n "${GTDBTK_DATA_PATH:-}" && -d "${GTDBTK_DATA_PATH}" ]] || die "GTDBTK_DATA_PATH is not set or is not a folder (GTDB-Tk reference package)."
    mkdir -p "${GTDBTK_DIR}"
    [[ -s "${GTDBTK_CACHE}" ]] || printf 'ASSEMBLY_ACCESSION\tGTDB_DOMAIN\tGTDB_PHYLUM\tGTDB_CLASS\tGTDB_ORDER\tGTDB_FAMILY\tGTDB_GENUS\tGTDB_SPECIES\tGTDB_TAXONOMY_SOURCE\n' > "${GTDBTK_CACHE}"

    # Genomes already classified by GTDB-Tk in a previous run are not run again
    awk 'NR == FNR { if (FNR > 1) done[$1] = 1; next } !($1 in done)' "${GTDBTK_CACHE}" "${WORK_DIR}/not_in_gtdb.txt" > "${WORK_DIR}/gtdbtk_pending.txt"

    : > "${WORK_DIR}/gtdbtk_batch.tsv"
    while IFS= read -r ACC; do
        FILE="$(genome_file_of "${CLEAN_DIR}" "${ACC}")"
        if [[ -n "${FILE}" ]]; then
            printf '%s\t%s\n' "$(absolute_path "${FILE}")" "${ACC}" >> "${WORK_DIR}/gtdbtk_batch.tsv"
        fi
    done < "${WORK_DIR}/gtdbtk_pending.txt"

    N_RUN="$(count_lines "${WORK_DIR}/gtdbtk_batch.tsv")"
    log "message" "GTDB-Tk: $(( N_NOT_IN_GTDB - $(count_lines "${WORK_DIR}/gtdbtk_pending.txt") )) genome(s) reused from cache; ${N_RUN} to classify."

    if [[ "${N_RUN}" -gt 0 ]]; then
        GTDBTK_VERSION="$(gtdbtk --version 2>&1 | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -n 1 || true)"
        RUN_DIR="${GTDBTK_DIR}/run_$(date +%Y%m%d_%H%M%S)"
        log "message" "Running GTDB-Tk ${GTDBTK_VERSION:-} classify_wf (${THREADS} CPUs). This step needs a large amount of RAM."
        # shellcheck disable=SC2086
        if ! gtdbtk classify_wf --batchfile "${WORK_DIR}/gtdbtk_batch.tsv" --out_dir "${RUN_DIR}" \
                --cpus "${THREADS}" --pplacer_cpus "$(( THREADS > 1 ? THREADS / 2 : 1 ))" \
                ${GTDBTK_EXTRA_ARGS:-} > "${RUN_DIR}.console.log" 2>&1; then
            log "error" "GTDB-Tk failed. Last lines of ${RUN_DIR}.console.log:"
            tail -n 20 "${RUN_DIR}.console.log" >&2
            exit 1
        fi
        parse_gtdbtk_summaries "${RUN_DIR}" "GTDB-Tk_${GTDBTK_VERSION:-unknown}" > "${WORK_DIR}/gtdbtk_new_rows.tsv"
        log "message" "GTDB-Tk classified $(count_lines "${WORK_DIR}/gtdbtk_new_rows.tsv") of ${N_RUN} genome(s)."

        # Genomes that GTDB-Tk could not classify are also cached (with NULL), so they are not run again
        awk -F'\t' -v OFS='\t' 'NR == FNR { done[$1] = 1; print; next }
            !($2 in done) { print $2, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "GTDB-Tk_UNCLASSIFIED" }' \
            "${WORK_DIR}/gtdbtk_new_rows.tsv" "${WORK_DIR}/gtdbtk_batch.tsv" > "${WORK_DIR}/gtdbtk_rows.tsv"
        { head -n 1 "${GTDBTK_CACHE}"; { cat "${WORK_DIR}/gtdbtk_rows.tsv"; tail -n +2 "${GTDBTK_CACHE}"; } | awk -F'\t' '!seen[$1]++'; } > "${GTDBTK_CACHE}.tmp"
        mv "${GTDBTK_CACHE}.tmp" "${GTDBTK_CACHE}"
    fi

    # Replace the rows of the genomes not found in GTDB with the GTDB-Tk result (when it is a classification)
    awk -F'\t' -v OFS='\t' 'NR == FNR { if (FNR > 1 && $2 != "NULL") tk[$1] = $0; next }
        FNR > 1 && $NF == "NULL" && ($1 in tk) { print tk[$1]; next } { print }' \
        "${GTDBTK_CACHE}" "${WORK_DIR}/taxonomy.tsv" > "${WORK_DIR}/taxonomy_final.tsv"
    mv "${WORK_DIR}/taxonomy_final.tsv" "${WORK_DIR}/taxonomy.tsv"
elif [[ "${N_NOT_IN_GTDB}" -gt 0 ]]; then
    log "message" "Genomes not found in the GTDB tables are listed in ${NON_CLASSIFIED_LIST} (use --gtdbtk to classify them with GTDB-Tk)."
fi

mv "${WORK_DIR}/taxonomy.tsv" "${TAXONOMY_TABLE}"
awk -F'\t' 'NR > 1 && $NF == "NULL" { print $1 }' "${TAXONOMY_TABLE}" > "${NON_CLASSIFIED_LIST}"

log "message" "Taxonomy saved to: ${TAXONOMY_TABLE}"
log "message" "Genomes without GTDB taxonomy (left out of the directory hierarchy): $(count_lines "${NON_CLASSIFIED_LIST}") -> ${NON_CLASSIFIED_LIST}"
