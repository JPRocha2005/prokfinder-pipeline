#!/usr/bin/env bash
#
# evaluate-genome-structure-seqkit.sh - Prokfinder module 4: structural evaluation of the genomes (seqkit)
#
# For each clean genome (small contigs already removed), seqkit stats measures the genome size, number of contigs,
# N50, number of ambiguous bases (N) and GC content. A genome passes (STRUCTURE_STATUS = PASSED) when:
#     number of contigs < --max-contigs (1000)    N50 >= --min-n50 (5000 bp)    ambiguous bases (N) < --max-n-bases (100000)
# Only these genomes go on to the gene prediction and the biological evaluation (CheckM2 and BUSCO).
#
# Re-execution: the seqkit results of each genome are kept in a cache table with the modification time of the genome
# file. The thresholds are applied on every run, so changing them does not run seqkit again.
#
# Usage: evaluate-genome-structure-seqkit.sh -l <accession_list> -i <clean_genomes_dir> -o <quality_dir> [-p <threads>]
#            [--max-contigs N] [--min-n50 BP] [--max-n-bases N]
#
# Output inside <quality_dir>:
#   genome-structure.tsv        ASSEMBLY_ACCESSION + seqkit metrics + STRUCTURE_STATUS + STRUCTURE_FAILED_CRITERIA
#   structure-passed-genomes.txt   accessions with STRUCTURE_STATUS = PASSED
#   seqkit-stats-cache.tsv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -i <clean_genomes_dir> -o <quality_dir> [-p <threads>] [thresholds]

  -l <file>             Accessions to evaluate (one per line)
  -i <dir>              Folder with the clean genomes (<ACCESSION>.fna or .fna.gz)
  -o <dir>              Output folder (tables and cache)
  -p <n>                CPUs used by seqkit (default: 4)
  --max-contigs <n>     The genome must have fewer contigs than this (default: 1000)
  --min-n50 <bp>        Minimum N50 (default: 5000)
  --max-n-bases <n>     The genome must have fewer ambiguous bases (N) than this (default: 100000)
EOF
    exit "${1:-1}"
}

ACC_LIST=""
CLEAN_DIR=""
QUALITY_DIR=""
THREADS=4
QC_MAX_CONTIGS=1000
QC_MIN_N50=5000
QC_MAX_N_BASES=100000

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -i) CLEAN_DIR="${2%/}"; shift 2 ;;
        -o) QUALITY_DIR="${2%/}"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        --max-contigs) QC_MAX_CONTIGS="$2"; shift 2 ;;
        --min-n50) QC_MIN_N50="$2"; shift 2 ;;
        --max-n-bases) QC_MAX_N_BASES="$2"; shift 2 ;;
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
for OPT in "--max-contigs:${QC_MAX_CONTIGS}" "--min-n50:${QC_MIN_N50}" "--max-n-bases:${QC_MAX_N_BASES}"; do
    [[ "${OPT#*:}" =~ ^[0-9]+$ ]] || die "${OPT%%:*} must be a non-negative integer (got '${OPT#*:}')."
done

mkdir -p "${QUALITY_DIR}"
setup_module_log "${QUALITY_DIR}/evaluate-genome-structure-seqkit.log"
require_programs seqkit

SEQKIT_CACHE="${QUALITY_DIR}/seqkit-stats-cache.tsv"
STRUCTURE_TABLE="${QUALITY_DIR}/genome-structure.tsv"
PASSED_LIST="${QUALITY_DIR}/structure-passed-genomes.txt"
[[ -s "${SEQKIT_CACHE}" ]] || printf 'ASSEMBLY_ACCESSION\tFILE_MTIME\tCLEAN_NUM_CONTIGS\tCLEAN_GENOME_SIZE_BP\tCLEAN_N50\tCLEAN_NUM_N_BASES\tCLEAN_GC_PERCENT\n' > "${SEQKIT_CACHE}"

make_temp_dir WORK_DIR
GENOMES="${WORK_DIR}/genomes.tsv"    # <accession> <TAB> <file> <TAB> <mtime>

# 1. Genome files of the requested accessions
while IFS= read -r ACC; do
    FILE="$(genome_file_of "${CLEAN_DIR}" "${ACC}")"
    if [[ -n "${FILE}" ]]; then
        printf '%s\t%s\t%s\n' "${ACC}" "$(absolute_path "${FILE}")" "$(file_mtime "${FILE}")"
    fi
done < <(tr -d '\r' < "${ACC_LIST}" | awk 'NF { print $1 }' | sort -u) > "${GENOMES}"

N_GENOMES="$(count_lines "${GENOMES}")"
log "message" "Structural evaluation (seqkit): ${N_GENOMES} clean genome(s) (${THREADS} CPUs)."
log "message" "  Criteria: contigs < ${QC_MAX_CONTIGS}; N50 >= ${QC_MIN_N50} bp; N bases < ${QC_MAX_N_BASES}"
[[ "${N_GENOMES}" -gt 0 ]] || die "No clean genome of the requested accessions was found in ${CLEAN_DIR}."

# 2. seqkit stats (only genomes not in the cache). "--gap-letters Nn" makes sum_gap count the ambiguous bases N/n.
pending_for_cache "${SEQKIT_CACHE}" "${GENOMES}" > "${WORK_DIR}/seqkit_pending.tsv"
N_PENDING="$(count_lines "${WORK_DIR}/seqkit_pending.tsv")"
log "message" "  $((N_GENOMES - N_PENDING)) genome(s) reused from cache; ${N_PENDING} to evaluate."

if [[ "${N_PENDING}" -gt 0 ]]; then
    log "message" "  seqkit: $(seqkit version 2>&1 | head -n 1)"
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

# 3. Result table: metrics of the cache + thresholds of this run
awk -F'\t' -v OFS='\t' -v cache="${SEQKIT_CACHE}" \
    -v max_contigs="${QC_MAX_CONTIGS}" -v min_n50="${QC_MIN_N50}" -v max_n="${QC_MAX_N_BASES}" '
    BEGIN {
        while ((getline line < cache) > 0) { split(line, f, "\t"); if (f[1] != "ASSEMBLY_ACCESSION") sk[f[1] SUBSEP f[2]] = line }
        print "ASSEMBLY_ACCESSION", "CLEAN_NUM_CONTIGS", "CLEAN_GENOME_SIZE_BP", "CLEAN_N50", "CLEAN_NUM_N_BASES", "CLEAN_GC_PERCENT", \
              "STRUCTURE_STATUS", "STRUCTURE_FAILED_CRITERIA"
    }
    {
        key = $1 SUBSEP $3
        if (!(key in sk)) { print $1, "NULL", "NULL", "NULL", "NULL", "NULL", "NOT_ASSESSED", "NULL"; next }
        split(sk[key], s, "\t")
        fails = ""
        if (s[3] + 0 >= max_contigs + 0) fails = fails ";CONTIGS>=" max_contigs
        if (s[5] + 0 < min_n50 + 0)      fails = fails ";N50<" min_n50
        if (s[6] + 0 >= max_n + 0)       fails = fails ";N_BASES>=" max_n
        if (fails != "") { sub(/^;/, "", fails); print $1, s[3], s[4], s[5], s[6], s[7], "FAILED", fails }
        else             print $1, s[3], s[4], s[5], s[6], s[7], "PASSED", "NULL"
    }
' "${GENOMES}" > "${STRUCTURE_TABLE}"

awk -F'\t' 'NR > 1 && $7 == "PASSED" { print $1 }' "${STRUCTURE_TABLE}" > "${PASSED_LIST}"
log "message" "Genomes that passed the structural evaluation: $(count_lines "${PASSED_LIST}") of ${N_GENOMES} (table: ${STRUCTURE_TABLE})"
