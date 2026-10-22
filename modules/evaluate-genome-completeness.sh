#!/usr/bin/env bash
#
# evaluate-genome-completeness.sh - Prokfinder module 6: biological evaluation of the genomes (CheckM2 and BUSCO)
#
# Both tools read the proteins predicted by predict-genes-pyrodigal.sh (no gene calling is repeated):
#   - CheckM2 (predict --genes): completeness and contamination. A genome passes (QC_STATUS = PASSED) when
#         completeness >= --min-completeness (75%)    contamination <= --max-contamination (5%)
#     (the genomes that failed the structural evaluation keep QC_STATUS = FAILED with its criteria)
#   - BUSCO (--mode proteins), only for the genomes that passed: single-copy orthologs complete (single or
#     duplicated), fragmented and missing. The lineage dataset follows the NCBI domain of the genome
#     (bacteria_odbNN or archaea_odbNN; odb12 with BUSCO >= 6, odb10 before). BUSCO is reported, not filtered.
#     The lineage datasets are downloaded once to the folder of -b (the reference databases folder).
#
# Re-execution: the CheckM2 and BUSCO results of each genome are kept in cache tables with the modification time of
# its protein file, so only new or changed genomes are evaluated again. The thresholds are applied on every run.
#
# Usage: evaluate-genome-completeness.sh -l <accession_list> -g <genes_dir> -s <genome-structure.tsv>
#            -d <genome-domain-ncbi.tsv> -o <quality_dir> -b <busco_download_dir> [-p <threads>]
#            [--min-completeness PCT] [--max-contamination PCT] [--skip-busco]
#
# Output inside <quality_dir>:
#   genome-quality.tsv         ASSEMBLY_ACCESSION CHECKM2_* BUSCO_* QC_STATUS QC_FAILED_CRITERIA (every genome of -s)
#   high-quality-genomes.txt   accessions with QC_STATUS = PASSED
#   checkm2-cache.tsv / busco-cache.tsv / checkm2-console.log / busco-console.log

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -g <genes_dir> -s <genome-structure.tsv> -d <genome-domain-ncbi.tsv>
       -o <quality_dir> -b <busco_download_dir> [-p <threads>] [thresholds] [--skip-busco]

  -l <file>                  Accessions to evaluate (genomes with predicted genes; one per line)
  -g <dir>                   Folder with the predicted proteins (<ACCESSION>.faa)
  -s <file>                  Structural evaluation table (genome-structure.tsv)
  -d <file>                  NCBI domain of each genome (genome-domain-ncbi.tsv), for the BUSCO lineage
  -o <dir>                   Output folder (tables and caches)
  -b <dir>                   Folder where BUSCO keeps its lineage datasets (downloaded on first use)
  -p <n>                     CPUs used by CheckM2 and BUSCO (default: 4)
  --min-completeness <pct>   Minimum CheckM2 completeness (default: 75)
  --max-contamination <pct>  Maximum CheckM2 contamination (default: 5)
  --skip-busco               Do not run BUSCO

CHECKM2DB must point to the CheckM2 database file (uniref100.KO.1.dmnd); the checkm2 module of the cluster sets it.
EOF
    exit "${1:-1}"
}

ACC_LIST=""
GENES_DIR=""
STRUCTURE_TABLE=""
DOMAIN_TABLE=""
QUALITY_DIR=""
BUSCO_DOWNLOADS=""
THREADS=4
QC_MIN_COMPLETENESS=75
QC_MAX_CONTAMINATION=5
SKIP_BUSCO=0
CHECKM2_BATCH_SIZE=1000   # genomes per CheckM2 run
BUSCO_BATCH_SIZE=500      # genomes per BUSCO run

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -g) GENES_DIR="${2%/}"; shift 2 ;;
        -s) STRUCTURE_TABLE="$2"; shift 2 ;;
        -d) DOMAIN_TABLE="$2"; shift 2 ;;
        -o) QUALITY_DIR="${2%/}"; shift 2 ;;
        -b) BUSCO_DOWNLOADS="${2%/}"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        --min-completeness) QC_MIN_COMPLETENESS="$2"; shift 2 ;;
        --max-contamination) QC_MAX_CONTAMINATION="$2"; shift 2 ;;
        --skip-busco) SKIP_BUSCO=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${ACC_LIST}" || -z "${GENES_DIR}" || -z "${STRUCTURE_TABLE}" || -z "${QUALITY_DIR}" ]]; then
    log "error" "Options -l, -g, -s and -o are required."
    usage
fi
if [[ "${SKIP_BUSCO}" -eq 0 && ( -z "${DOMAIN_TABLE}" || -z "${BUSCO_DOWNLOADS}" ) ]]; then
    die "Options -d and -b are required for BUSCO (or use --skip-busco)."
fi
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"
[[ -f "${STRUCTURE_TABLE}" ]] || die "Structural evaluation table not found: ${STRUCTURE_TABLE}"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "The number of threads must be a positive integer (got '${THREADS}')."
for OPT in "--min-completeness:${QC_MIN_COMPLETENESS}" "--max-contamination:${QC_MAX_CONTAMINATION}"; do
    if [[ ! "${OPT#*:}" =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v v="${OPT#*:}" 'BEGIN { exit !(v + 0 <= 100) }'; then
        die "${OPT%%:*} must be a percentage between 0 and 100 (got '${OPT#*:}')."
    fi
done

mkdir -p "${QUALITY_DIR}"
setup_module_log "${QUALITY_DIR}/evaluate-genome-completeness.log"
require_programs checkm2
[[ "${SKIP_BUSCO}" -eq 1 ]] || require_programs busco
if [[ -z "${CHECKM2DB:-}" ]]; then
    die "CHECKM2DB is not set. It must be the path of the CheckM2 database file (uniref100.KO.1.dmnd)."
elif [[ ! -f "${CHECKM2DB}" ]]; then
    die "CheckM2 database not found: ${CHECKM2DB}"
fi

CHECKM2_CACHE="${QUALITY_DIR}/checkm2-cache.tsv"
BUSCO_CACHE="${QUALITY_DIR}/busco-cache.tsv"
CHECKM2_CONSOLE_LOG="${QUALITY_DIR}/checkm2-console.log"
BUSCO_CONSOLE_LOG="${QUALITY_DIR}/busco-console.log"
QUALITY_TABLE="${QUALITY_DIR}/genome-quality.tsv"
HIGH_QUALITY_LIST="${QUALITY_DIR}/high-quality-genomes.txt"
# (the cache of the version that ran CheckM2 on the genomes is not compatible: CHECKM2_INPUT=proteins marks this one)
[[ -s "${CHECKM2_CACHE}" ]] && head -n 1 "${CHECKM2_CACHE}" | grep -q "CHECKM2_INPUT" || \
    printf 'ASSEMBLY_ACCESSION\tFILE_MTIME\tCHECKM2_INPUT\tCHECKM2_COMPLETENESS\tCHECKM2_CONTAMINATION\tCHECKM2_MODEL\n' > "${CHECKM2_CACHE}"
[[ -s "${BUSCO_CACHE}" ]] || printf 'ASSEMBLY_ACCESSION\tFILE_MTIME\tBUSCO_LINEAGE\tBUSCO_COMPLETENESS\tBUSCO_COMPLETE\tBUSCO_SINGLE\tBUSCO_DUPLICATED\tBUSCO_FRAGMENTED\tBUSCO_MISSING\tBUSCO_TOTAL\n' > "${BUSCO_CACHE}"

make_temp_dir WORK_DIR
PROTEINS="${WORK_DIR}/proteins.tsv"    # <accession> <TAB> <.faa> <TAB> <mtime>

# 1. Protein files of the requested accessions
while IFS= read -r ACC; do
    FILE="${GENES_DIR}/${ACC}.faa"
    if [[ -s "${FILE}" ]]; then
        printf '%s\t%s\t%s\n' "${ACC}" "$(absolute_path "${FILE}")" "$(file_mtime "${FILE}")"
    fi
done < <(tr -d '\r' < "${ACC_LIST}" | awk 'NF { print $1 }' | sort -u) > "${PROTEINS}"
N_GENOMES="$(count_lines "${PROTEINS}")"
log "message" "Biological evaluation: ${N_GENOMES} genome(s) with predicted proteins (${THREADS} CPUs)."
log "message" "  Criteria: CheckM2 completeness >= ${QC_MIN_COMPLETENESS}%; contamination <= ${QC_MAX_CONTAMINATION}%$([[ "${SKIP_BUSCO}" -eq 0 ]] && echo "; BUSCO is reported (not filtered)")"

# Links the files of <rows> (<accession> <TAB> <file> ...) as <dir>/<accession>.<ext>
link_inputs() {
    local ROWS="$1" DIR="$2" EXT="$3" ACC FILE REST
    mkdir -p "${DIR}"
    while IFS=$'\t' read -r ACC FILE REST; do
        ln -sf "${FILE}" "${DIR}/${ACC}.${EXT}"
    done < "${ROWS}"
}

# 2. CheckM2 on the predicted proteins (--genes), genomes not in the cache
pending_for_cache "${CHECKM2_CACHE}" "${PROTEINS}" > "${WORK_DIR}/checkm2_pending.tsv"
N_PENDING="$(count_lines "${WORK_DIR}/checkm2_pending.tsv")"
log "message" "CheckM2 (--genes): $((N_GENOMES - N_PENDING)) genome(s) reused from cache; ${N_PENDING} to evaluate."

if [[ "${N_PENDING}" -gt 0 ]]; then
    log "message" "  checkm2: $(run_in_own_env checkm2 --version 2>&1 | head -n 1); database: ${CHECKM2DB}"
    split -l "${CHECKM2_BATCH_SIZE}" -d -a 4 "${WORK_DIR}/checkm2_pending.tsv" "${WORK_DIR}/checkm2_batch_"
    for BATCH in "${WORK_DIR}"/checkm2_batch_*; do
        BATCH_NAME="$(basename "${BATCH}")"
        log "message" "  CheckM2 ${BATCH_NAME#checkm2_}: $(count_lines "${BATCH}") genome(s)..."
        INPUT_DIR="${WORK_DIR}/${BATCH_NAME}_input"
        OUT_DIR="${WORK_DIR}/${BATCH_NAME}_out"
        link_inputs "${BATCH}" "${INPUT_DIR}" faa
        if ! run_in_own_env checkm2 predict --genes --input "${INPUT_DIR}" --extension faa --output-directory "${OUT_DIR}" \
                --threads "${THREADS}" --database_path "${CHECKM2DB}" --remove_intermediates --force \
                >> "${CHECKM2_CONSOLE_LOG}" 2>&1; then
            log "error" "checkm2 predict failed. Last lines of ${CHECKM2_CONSOLE_LOG}:"
            tail -n 20 "${CHECKM2_CONSOLE_LOG}" >&2
            exit 1
        fi
        [[ -s "${OUT_DIR}/quality_report.tsv" ]] || die "CheckM2 did not write quality_report.tsv. See ${CHECKM2_CONSOLE_LOG}."
        awk -F'\t' -v OFS='\t' -v batch="${BATCH}" '
            BEGIN { while ((getline line < batch) > 0) { split(line, f, "\t"); mtime_of[f[1]] = f[3] } }
            NR == 1 {
                for (i = 1; i <= NF; i++) col[$i] = i
                if (!("Name" in col) || !("Completeness" in col) || !("Contamination" in col)) {
                    print "checkm2 quality_report.tsv: the columns Name, Completeness and Contamination are required" > "/dev/stderr"
                    exit 2
                }
                next
            }
            {
                if (!match($col["Name"], /GC[AF]_[0-9]+[.][0-9]+/)) next
                acc = substr($col["Name"], RSTART, RLENGTH)
                model = ("Completeness_Model_Used" in col) ? $col["Completeness_Model_Used"] : "NULL"
                print acc, mtime_of[acc], "proteins", $col["Completeness"], $col["Contamination"], model
            }
        ' "${OUT_DIR}/quality_report.tsv" > "${WORK_DIR}/checkm2_new_rows.tsv" \
            || die "Could not read the CheckM2 report (its columns may have changed)."
        # The cache is updated after each batch: an interrupted run keeps the batches already evaluated
        update_cache "${CHECKM2_CACHE}" "${WORK_DIR}/checkm2_new_rows.tsv"
        rm -rf "${INPUT_DIR}" "${OUT_DIR}"
    done
fi

# Genomes that pass the CheckM2 criteria (with the result of their current protein file)
awk -F'\t' -v min_comp="${QC_MIN_COMPLETENESS}" -v max_cont="${QC_MAX_CONTAMINATION}" '
    NR == FNR { if (FNR > 1) { comp[$1 SUBSEP $2] = $4; cont[$1 SUBSEP $2] = $5 }; next }
    (($1 SUBSEP $3) in comp) && comp[$1 SUBSEP $3] + 0 >= min_comp + 0 && cont[$1 SUBSEP $3] + 0 <= max_cont + 0
' "${CHECKM2_CACHE}" "${PROTEINS}" > "${WORK_DIR}/checkm2_passed.tsv"

# 3. BUSCO (proteins mode) on the genomes that passed CheckM2, one run per NCBI domain (lineage)
if [[ "${SKIP_BUSCO}" -eq 0 && -s "${WORK_DIR}/checkm2_passed.tsv" ]]; then
    BUSCO_VERSION="$(run_in_own_env busco --version 2>&1 | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -n 1 || true)"
    ODB="odb10"
    [[ "${BUSCO_VERSION%%.*}" =~ ^[0-9]+$ ]] && (( ${BUSCO_VERSION%%.*} >= 6 )) && ODB="odb12"

    # <accession> <file> <mtime> <lineage> (genomes of unknown domain use the bacteria dataset)
    awk -F'\t' -v OFS='\t' -v odb="${ODB}" '
        NR == FNR { if ($0 !~ /^#/ && $1 != "ASSEMBLY_ACCESSION") dom[$1] = $2; next }
        { print $0, ((dom[$1] == "ARCHAEA") ? "archaea_" : "bacteria_") odb }
    ' "${DOMAIN_TABLE}" "${WORK_DIR}/checkm2_passed.tsv" > "${WORK_DIR}/busco_candidates.tsv"
    # pending: not in the cache with the same mtime and lineage
    awk -F'\t' 'NR == FNR { if (FNR > 1) done[$1 SUBSEP $2 SUBSEP $3] = 1; next } !(($1 SUBSEP $3 SUBSEP $4) in done)' \
        "${BUSCO_CACHE}" "${WORK_DIR}/busco_candidates.tsv" > "${WORK_DIR}/busco_pending.tsv"
    N_CANDIDATES="$(count_lines "${WORK_DIR}/busco_candidates.tsv")"
    N_PENDING="$(count_lines "${WORK_DIR}/busco_pending.tsv")"
    log "message" "BUSCO ${BUSCO_VERSION:-} (proteins): ${N_CANDIDATES} genome(s) passed CheckM2; $((N_CANDIDATES - N_PENDING)) reused from cache; ${N_PENDING} to evaluate."

    mkdir -p "${BUSCO_DOWNLOADS}"
    for LINEAGE in $(cut -f4 "${WORK_DIR}/busco_pending.tsv" | sort -u); do
        awk -F'\t' -v l="${LINEAGE}" '$4 == l' "${WORK_DIR}/busco_pending.tsv" > "${WORK_DIR}/busco_${LINEAGE}.tsv"
        split -l "${BUSCO_BATCH_SIZE}" -d -a 4 "${WORK_DIR}/busco_${LINEAGE}.tsv" "${WORK_DIR}/busco_${LINEAGE}_batch_"
        for BATCH in "${WORK_DIR}/busco_${LINEAGE}_batch_"*; do
            BATCH_NAME="$(basename "${BATCH}")"
            log "message" "  BUSCO ${LINEAGE}, ${BATCH_NAME##*_batch_}: $(count_lines "${BATCH}") genome(s)..."
            INPUT_DIR="${WORK_DIR}/${BATCH_NAME}_input"
            RUN_NAME="${BATCH_NAME}_run"
            link_inputs "${BATCH}" "${INPUT_DIR}" faa
            # (BUSCO runs in a temporary folder: it writes relative paths and many files per genome)
            if ! ( cd "${WORK_DIR}" && run_in_own_env busco --in "${INPUT_DIR}" --mode proteins --lineage_dataset "${LINEAGE}" \
                    --out "${RUN_NAME}" --out_path "${WORK_DIR}" --download_path "${BUSCO_DOWNLOADS}" \
                    --cpu "${THREADS}" --force --quiet ) >> "${BUSCO_CONSOLE_LOG}" 2>&1; then
                log "warning" "BUSCO failed for ${BATCH_NAME} (see ${BUSCO_CONSOLE_LOG}); its genomes have no BUSCO result."
                continue
            fi
            # short_summary.specific.<lineage>.<input>.txt of each genome: counts of complete (single, duplicated),
            # fragmented and missing orthologs
            find "${WORK_DIR}/${RUN_NAME}" -type f -name 'short_summary.specific.*.txt' -print0 \
                | xargs -0 -r awk -v OFS='\t' -v batch="${BATCH}" '
                    BEGIN { while ((getline line < batch) > 0) { split(line, f, "\t"); mtime_of[f[1]] = f[3]; lin[f[1]] = f[4] } }
                    FNR == 1 {
                        if (acc != "") print acc, mtime_of[acc], lin[acc], pct, c, s, d, fr, m, n
                        acc = ""; pct = c = s = d = fr = m = n = "NULL"
                        if (match(FILENAME, /GC[AF]_[0-9]+[.][0-9]+/)) acc = substr(FILENAME, RSTART, RLENGTH)
                    }
                    /C:[0-9.]+%/ { if (match($0, /C:[0-9.]+/)) pct = substr($0, RSTART + 2, RLENGTH - 2) }
                    /Complete BUSCOs/                    { c = $1 }
                    /Complete and single-copy BUSCOs/    { s = $1 }
                    /Complete and duplicated BUSCOs/     { d = $1 }
                    /Fragmented BUSCOs/                  { fr = $1 }
                    /Missing BUSCOs/                     { m = $1 }
                    /Total BUSCO groups searched/        { n = $1 }
                    END { if (acc != "") print acc, mtime_of[acc], lin[acc], pct, c, s, d, fr, m, n }
                ' > "${WORK_DIR}/busco_new_rows.tsv"
            update_cache "${BUSCO_CACHE}" "${WORK_DIR}/busco_new_rows.tsv"
            rm -rf "${INPUT_DIR}" "${WORK_DIR:?}/${RUN_NAME}"
        done
    done
fi

# 4. Final table: every genome of the structural evaluation, with the results of the caches and this run's thresholds
awk -F'\t' -v OFS='\t' -v proteins="${PROTEINS}" -v checkm2_cache="${CHECKM2_CACHE}" -v busco_cache="${BUSCO_CACHE}" \
    -v min_comp="${QC_MIN_COMPLETENESS}" -v max_cont="${QC_MAX_CONTAMINATION}" -v skip_busco="${SKIP_BUSCO}" '
    BEGIN {
        while ((getline line < proteins) > 0) { split(line, f, "\t"); mtime[f[1]] = f[3] }
        while ((getline line < checkm2_cache) > 0) { split(line, f, "\t"); if (f[1] != "ASSEMBLY_ACCESSION") cm[f[1] SUBSEP f[2]] = f[4] OFS f[5] OFS f[6] }
        while ((getline line < busco_cache) > 0) {
            split(line, f, "\t")
            if (f[1] != "ASSEMBLY_ACCESSION") { sub(/^[^\t]*\t[^\t]*\t/, "", line); bu[f[1] SUBSEP f[2]] = line }
        }
        print "ASSEMBLY_ACCESSION", "CHECKM2_COMPLETENESS", "CHECKM2_CONTAMINATION", "CHECKM2_MODEL", "BUSCO_LINEAGE",
              "BUSCO_COMPLETENESS", "BUSCO_COMPLETE", "BUSCO_SINGLE", "BUSCO_DUPLICATED", "BUSCO_FRAGMENTED", "BUSCO_MISSING",
              "BUSCO_TOTAL", "QC_STATUS", "QC_FAILED_CRITERIA"
    }
    FNR == 1 { for (i = 1; i <= NF; i++) col[$i] = i; next }
    {
        acc = $1
        checkm2 = "NULL" OFS "NULL" OFS "NULL"
        busco = "NULL" OFS "NULL" OFS "NULL" OFS "NULL" OFS "NULL" OFS "NULL" OFS "NULL" OFS "NULL"
        if ($col["STRUCTURE_STATUS"] != "PASSED") {
            status = ($col["STRUCTURE_STATUS"] == "FAILED") ? "FAILED" : $col["STRUCTURE_STATUS"]
            print acc, checkm2, busco, status, $col["STRUCTURE_FAILED_CRITERIA"]
            next
        }
        if (!(acc in mtime)) { print acc, checkm2, busco, "NO_PREDICTED_GENES", "NULL"; next }
        key = acc SUBSEP mtime[acc]
        if (key in bu) busco = bu[key]
        if (!(key in cm)) { print acc, checkm2, busco, "NO_CHECKM2_RESULT", "NULL"; next }
        checkm2 = cm[key]
        split(checkm2, c, OFS)
        fails = ""
        if (c[1] + 0 < min_comp + 0) fails = fails ";COMPLETENESS<" min_comp
        if (c[2] + 0 > max_cont + 0) fails = fails ";CONTAMINATION>" max_cont
        if (fails != "") { sub(/^;/, "", fails); print acc, checkm2, busco, "FAILED", fails }
        else print acc, checkm2, busco, "PASSED", "NULL"
    }
' "${STRUCTURE_TABLE}" > "${QUALITY_TABLE}"

awk -F'\t' 'NR > 1 && $13 == "PASSED" { print $1 }' "${QUALITY_TABLE}" > "${HIGH_QUALITY_LIST}"
log "message" "Genome quality overview (genomes per result):"
while IFS=$'\t' read -r STATUS COUNT; do
    log "message" "    ${STATUS}: ${COUNT}"
done < <(awk -F'\t' 'NR > 1 { c[$13]++ } END { for (s in c) printf "%s\t%d\n", s, c[s] }' "${QUALITY_TABLE}" | LC_ALL=C sort)
log "message" "High-quality genomes: $(count_lines "${HIGH_QUALITY_LIST}") (list: ${HIGH_QUALITY_LIST}; table: ${QUALITY_TABLE})"
