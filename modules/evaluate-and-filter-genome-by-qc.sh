#!/usr/bin/env bash
#
# evaluate-and-filter-genome-by-qc.sh - Prokfinder module 4: genome quality assessment and filtering
#
# For each clean genome (contigs < 1000 bp already removed):
#   - seqkit stats: genome size, number of contigs, N50, number of ambiguous bases (N) and GC content
#   - CheckM2: completeness and contamination (only for the genomes that pass the seqkit criteria)
# A genome is kept (QC_STATUS = PASSED) when:
#     number of contigs < QC_MAX_CONTIGS (1000)        N50 >= QC_MIN_N50 (5000 bp)
#     ambiguous bases (N) < QC_MAX_N_BASES (100000)    CheckM2 completeness >= QC_MIN_COMPLETENESS (75%)
#     CheckM2 contamination <= QC_MAX_CONTAMINATION (5%)
# The thresholds can be changed through these environment variables.
#
# Re-execution: the seqkit and CheckM2 results of each genome are kept in cache tables together with the
# modification time of the genome file, so only new or changed genomes are assessed again. The thresholds are
# applied on every run, so changing them does not require re-running seqkit or CheckM2.
#
# Usage: evaluate-and-filter-genome-by-qc.sh -l <accession_list> -i <clean_genomes_dir> -o <quality_dir> [-p <threads>]
#
# Output inside <quality_dir>:
#   genome-quality.tsv          ASSEMBLY_ACCESSION + metrics + QC_STATUS + QC_FAILED_CRITERIA
#   high-quality-genomes.txt    accessions with QC_STATUS = PASSED
#   seqkit-stats-cache.tsv / checkm2-cache.tsv / checkm2-console.log

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -i <clean_genomes_dir> -o <quality_dir> [-p <threads>]

  -l <file>   Accessions to assess (one per line)
  -i <dir>    Folder with the clean genomes (<ACCESSION>.fna or .fna.gz)
  -o <dir>    Output folder (tables and caches)
  -p <n>      CPUs used by seqkit and CheckM2 (default: number of CPUs)

Environment variables (thresholds): QC_MAX_CONTIGS (1000), QC_MIN_N50 (5000), QC_MAX_N_BASES (100000),
QC_MIN_COMPLETENESS (75), QC_MAX_CONTAMINATION (5), CHECKM2_BATCH_SIZE (genomes per CheckM2 run, 500).
CHECKM2DB must point to the CheckM2 database file (uniref100.KO.1.dmnd).
EOF
    exit "${1:-1}"
}

ACC_LIST=""
CLEAN_DIR=""
QUALITY_DIR=""
THREADS="$(nproc 2>/dev/null || echo 4)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -i) CLEAN_DIR="${2%/}"; shift 2 ;;
        -o) QUALITY_DIR="${2%/}"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${ACC_LIST}" || -z "${CLEAN_DIR}" || -z "${QUALITY_DIR}" ]]; then
    log "error" "Options -l, -i and -o are required."
    usage
fi
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "The number of threads must be a positive integer (got '${THREADS}')."

QC_MAX_CONTIGS="${QC_MAX_CONTIGS:-1000}"
QC_MIN_N50="${QC_MIN_N50:-5000}"
QC_MAX_N_BASES="${QC_MAX_N_BASES:-100000}"
QC_MIN_COMPLETENESS="${QC_MIN_COMPLETENESS:-75}"
QC_MAX_CONTAMINATION="${QC_MAX_CONTAMINATION:-5}"
CHECKM2_BATCH_SIZE="${CHECKM2_BATCH_SIZE:-500}"

for VAR in QC_MAX_CONTIGS QC_MIN_N50 QC_MAX_N_BASES CHECKM2_BATCH_SIZE; do
    [[ "${!VAR}" =~ ^[0-9]+$ ]] || die "${VAR} must be a non-negative integer (got '${!VAR}')."
done
for VAR in QC_MIN_COMPLETENESS QC_MAX_CONTAMINATION; do
    if [[ ! "${!VAR}" =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v v="${!VAR}" 'BEGIN { exit !(v + 0 <= 100) }'; then
        die "${VAR} must be a percentage between 0 and 100 (got '${!VAR}')."
    fi
done

mkdir -p "${QUALITY_DIR}"
setup_module_log "${QUALITY_DIR}/evaluate-and-filter-genome-by-qc.log"
require_programs seqkit checkm2

if [[ -z "${CHECKM2DB:-}" ]]; then
    die "CHECKM2DB is not set. It must be the path of the CheckM2 database file (uniref100.KO.1.dmnd)."
elif [[ ! -f "${CHECKM2DB}" ]]; then
    die "CheckM2 database not found: ${CHECKM2DB}"
fi

SEQKIT_CACHE="${QUALITY_DIR}/seqkit-stats-cache.tsv"
CHECKM2_CACHE="${QUALITY_DIR}/checkm2-cache.tsv"
CHECKM2_CONSOLE_LOG="${QUALITY_DIR}/checkm2-console.log"
QUALITY_TABLE="${QUALITY_DIR}/genome-quality.tsv"
HIGH_QUALITY_LIST="${QUALITY_DIR}/high-quality-genomes.txt"

[[ -s "${SEQKIT_CACHE}" ]] || printf 'ASSEMBLY_ACCESSION\tFILE_MTIME\tCLEAN_NUM_CONTIGS\tCLEAN_GENOME_SIZE_BP\tCLEAN_N50\tCLEAN_NUM_N_BASES\tCLEAN_GC_PERCENT\n' > "${SEQKIT_CACHE}"
[[ -s "${CHECKM2_CACHE}" ]] || printf 'ASSEMBLY_ACCESSION\tFILE_MTIME\tCHECKM2_COMPLETENESS\tCHECKM2_CONTAMINATION\tCHECKM2_MODEL\tCHECKM2_CODING_DENSITY\n' > "${CHECKM2_CACHE}"

make_temp_dir WORK_DIR
GENOMES="${WORK_DIR}/genomes.tsv"    # <accession> <TAB> <file> <TAB> <mtime>

# Prints the rows of <genomes.tsv> whose (accession, mtime) is not in <cache>
pending_for_cache() {
    awk -F'\t' 'NR == FNR { if (FNR > 1) done[$1 SUBSEP $2] = 1; next } !(($1 SUBSEP $3) in done)' "$1" "$2"
}

# Appends <new_rows> to <cache>, replacing the old rows of the same accessions
update_cache() {
    local CACHE="$1"
    local NEW_ROWS="$2"
    { head -n 1 "${CACHE}"; { cat "${NEW_ROWS}"; tail -n +2 "${CACHE}"; } | awk -F'\t' '!seen[$1]++' | LC_ALL=C sort -t $'\t' -k1,1; } > "${CACHE}.tmp"
    mv "${CACHE}.tmp" "${CACHE}"
}

# 1. Genome files of the requested accessions
while IFS= read -r ACC; do
    FILE="$(genome_file_of "${CLEAN_DIR}" "${ACC}")"
    if [[ -n "${FILE}" ]]; then
        printf '%s\t%s\t%s\n' "${ACC}" "$(absolute_path "${FILE}")" "$(file_mtime "${FILE}")"
    fi
done < <(tr -d '\r' < "${ACC_LIST}" | awk 'NF { print $1 }' | sort -u) > "${GENOMES}"

N_GENOMES="$(count_lines "${GENOMES}")"
log "message" "Quality assessment: ${N_GENOMES} clean genome(s) (${THREADS} CPUs)."
log "message" "  Criteria: contigs < ${QC_MAX_CONTIGS}; N50 >= ${QC_MIN_N50} bp; N bases < ${QC_MAX_N_BASES}; completeness >= ${QC_MIN_COMPLETENESS}%; contamination <= ${QC_MAX_CONTAMINATION}%"
if [[ "${N_GENOMES}" -eq 0 ]]; then
    die "No clean genome of the requested accessions was found in ${CLEAN_DIR}."
fi

# 2. seqkit stats (only genomes not in the cache). "-G Nn" makes sum_gap count the ambiguous bases N/n.
pending_for_cache "${SEQKIT_CACHE}" "${GENOMES}" > "${WORK_DIR}/seqkit_pending.tsv"
N_PENDING="$(count_lines "${WORK_DIR}/seqkit_pending.tsv")"
log "message" "Step 1 (seqkit stats): $((N_GENOMES - N_PENDING)) genome(s) reused from cache; ${N_PENDING} to assess."

if [[ "${N_PENDING}" -gt 0 ]]; then
    log "message" "seqkit: $(seqkit version 2>&1 | head -n 1)"
    cut -f2 "${WORK_DIR}/seqkit_pending.tsv" > "${WORK_DIR}/seqkit_files.txt"
    seqkit stats --quiet --tabular --all --gap-letters "Nn" --threads "${THREADS}" \
        --infile-list "${WORK_DIR}/seqkit_files.txt" > "${WORK_DIR}/seqkit_stats.tsv" \
        || die "seqkit stats failed."

    awk -F'\t' -v OFS='\t' -v pending="${WORK_DIR}/seqkit_pending.tsv" '
        BEGIN { while ((getline line < pending) > 0) { split(line, f, "\t"); acc_of[f[2]] = f[1]; mtime_of[f[2]] = f[3] } }
        NR == 1 {
            for (i = 1; i <= NF; i++) col[$i] = i
            if (!("file" in col) || !("num_seqs" in col) || !("sum_len" in col) || !("N50" in col) || !("sum_gap" in col) || !("GC(%)" in col)) {
                print "seqkit stats: the columns file, num_seqs, sum_len, N50, sum_gap and GC(%) are required" > "/dev/stderr"
                exit 2
            }
            next
        }
        ($col["file"] in acc_of) {
            fn = $col["file"]
            print acc_of[fn], mtime_of[fn], $col["num_seqs"], $col["sum_len"], $col["N50"], $col["sum_gap"], $col["GC(%)"]
        }
    ' "${WORK_DIR}/seqkit_stats.tsv" > "${WORK_DIR}/seqkit_new_rows.tsv" || die "Could not read the seqkit stats table."
    update_cache "${SEQKIT_CACHE}" "${WORK_DIR}/seqkit_new_rows.tsv"
fi

# 3. CheckM2 only for the genomes that pass the seqkit criteria and are not in the cache
awk -F'\t' -v max_contigs="${QC_MAX_CONTIGS}" -v min_n50="${QC_MIN_N50}" -v max_n="${QC_MAX_N_BASES}" '
    NR == FNR { if (FNR > 1) { mt[$1] = $2; ok[$1] = ($3 + 0 < max_contigs + 0 && $5 + 0 >= min_n50 + 0 && $6 + 0 < max_n + 0) }; next }
    ($1 in ok) && ok[$1] && mt[$1] == $3
' "${SEQKIT_CACHE}" "${GENOMES}" > "${WORK_DIR}/checkm2_candidates.tsv"

pending_for_cache "${CHECKM2_CACHE}" "${WORK_DIR}/checkm2_candidates.tsv" > "${WORK_DIR}/checkm2_pending.tsv"
N_CANDIDATES="$(count_lines "${WORK_DIR}/checkm2_candidates.tsv")"
N_PENDING="$(count_lines "${WORK_DIR}/checkm2_pending.tsv")"
log "message" "Step 2 (CheckM2): ${N_CANDIDATES} genome(s) passed the seqkit criteria; $((N_CANDIDATES - N_PENDING)) reused from cache; ${N_PENDING} to assess."

if [[ "${N_PENDING}" -gt 0 ]]; then
    log "message" "checkm2: $(checkm2 --version 2>&1 | head -n 1); database: ${CHECKM2DB}"
    split -l "${CHECKM2_BATCH_SIZE}" -d -a 4 "${WORK_DIR}/checkm2_pending.tsv" "${WORK_DIR}/checkm2_batch_"

    for BATCH in "${WORK_DIR}"/checkm2_batch_*; do
        BATCH_NAME="$(basename "${BATCH}")"
        log "message" "  CheckM2 ${BATCH_NAME#checkm2_}: $(count_lines "${BATCH}") genome(s)..."
        : > "${WORK_DIR}/checkm2_new_rows.tsv"

        # CheckM2 accepts a single extension per run: plain and compressed genomes are run separately
        for EXT in fna fna.gz; do
            INPUT_DIR="${WORK_DIR}/run_${BATCH_NAME}_${EXT//./_}_input"
            OUT_DIR="${WORK_DIR}/run_${BATCH_NAME}_${EXT//./_}_out"
            mkdir -p "${INPUT_DIR}"
            awk -F'\t' -v gz="$([[ "${EXT}" == "fna.gz" ]] && echo 1 || echo 0)" '($2 ~ /[.]gz$/) == (gz + 0)' "${BATCH}" \
                | while IFS=$'\t' read -r ACC FILE _; do ln -sf "${FILE}" "${INPUT_DIR}/${ACC}.${EXT}"; done
            if [[ -z "$(ls -A "${INPUT_DIR}")" ]]; then
                continue
            fi

            if ! checkm2 predict --input "${INPUT_DIR}" --output-directory "${OUT_DIR}" --extension "${EXT}" \
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
                    density = ("Coding_Density" in col) ? $col["Coding_Density"] : "NULL"
                    print acc, mtime_of[acc], $col["Completeness"], $col["Contamination"], model, density
                }
            ' "${OUT_DIR}/quality_report.tsv" >> "${WORK_DIR}/checkm2_new_rows.tsv" \
                || die "Could not read the CheckM2 report (its columns may have changed)."
        done

        # The cache is updated after each batch: an interrupted run keeps the batches already assessed
        update_cache "${CHECKM2_CACHE}" "${WORK_DIR}/checkm2_new_rows.tsv"
    done
fi

# 4. Final table: metrics of the caches + thresholds of this run
awk -F'\t' -v OFS='\t' -v seqkit_cache="${SEQKIT_CACHE}" -v checkm2_cache="${CHECKM2_CACHE}" \
    -v max_contigs="${QC_MAX_CONTIGS}" -v min_n50="${QC_MIN_N50}" -v max_n="${QC_MAX_N_BASES}" \
    -v min_comp="${QC_MIN_COMPLETENESS}" -v max_cont="${QC_MAX_CONTAMINATION}" '
    BEGIN {
        while ((getline line < seqkit_cache) > 0) { split(line, f, "\t"); if (f[1] != "ASSEMBLY_ACCESSION") sk[f[1] SUBSEP f[2]] = line }
        while ((getline line < checkm2_cache) > 0) { split(line, f, "\t"); if (f[1] != "ASSEMBLY_ACCESSION") cm[f[1] SUBSEP f[2]] = line }
        print "ASSEMBLY_ACCESSION", "CLEAN_NUM_CONTIGS", "CLEAN_GENOME_SIZE_BP", "CLEAN_N50", "CLEAN_NUM_N_BASES", "CLEAN_GC_PERCENT", \
              "CHECKM2_COMPLETENESS", "CHECKM2_CONTAMINATION", "CHECKM2_MODEL", "QC_STATUS", "QC_FAILED_CRITERIA"
    }
    {
        acc = $1; key = $1 SUBSEP $3
        if (!(key in sk)) { print acc, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "NOT_ASSESSED", "NULL"; next }
        split(sk[key], s, "\t")
        fails = ""
        if (s[3] + 0 >= max_contigs + 0) fails = fails ";CONTIGS>=" max_contigs
        if (s[5] + 0 < min_n50 + 0)      fails = fails ";N50<" min_n50
        if (s[6] + 0 >= max_n + 0)       fails = fails ";N_BASES>=" max_n
        comp = "NULL"; cont = "NULL"; model = "NULL"
        if (key in cm) {
            split(cm[key], c, "\t"); comp = c[3]; cont = c[4]; model = c[5]
            if (comp + 0 < min_comp + 0) fails = fails ";COMPLETENESS<" min_comp
            if (cont + 0 > max_cont + 0) fails = fails ";CONTAMINATION>" max_cont
        }
        if (fails != "")        { status = "FAILED"; sub(/^;/, "", fails) }
        else if (!(key in cm))  { status = "NO_CHECKM2_RESULT"; fails = "NULL" }
        else                    { status = "PASSED"; fails = "NULL" }
        print acc, s[3], s[4], s[5], s[6], s[7], comp, cont, model, status, fails
    }
' "${GENOMES}" > "${QUALITY_TABLE}"

awk -F'\t' 'NR > 1 && $10 == "PASSED" { print $1 }' "${QUALITY_TABLE}" > "${HIGH_QUALITY_LIST}"

log "message" "Quality assessment overview (genomes per result):"
while IFS=$'\t' read -r STATUS COUNT; do
    log "message" "    ${STATUS}: ${COUNT}"
done < <(awk -F'\t' 'NR > 1 { c[$10]++ } END { for (s in c) printf "%s\t%d\n", s, c[s] }' "${QUALITY_TABLE}" | LC_ALL=C sort)
log "message" "High-quality genomes: $(count_lines "${HIGH_QUALITY_LIST}") of ${N_GENOMES} (list: ${HIGH_QUALITY_LIST}; table: ${QUALITY_TABLE})"
