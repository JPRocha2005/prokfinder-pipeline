#!/usr/bin/env bash
#
# dereplicate-genomes.sh - Prokfinder module 9: dereplication of the genomes of each taxon (skDER or CiDDER)
#
# The dereplication runs inside each leaf directory of the taxonomic hierarchy, i.e. only genomes with the same
# lineage (same species, or the lowest known rank) are compared. Leaves with a single genome need no dereplication:
# the genome is its own representative. Depending on the parameters, a leaf may keep more than one representative
# genome when there is divergence/variability inside the taxon. Two methods (both from the skDER package):
#   skder  (default) ANI between the genomes (skani): genomes above the ANI cutoff (99.5% by default) with an aligned
#          fraction of 90% are redundant; --mode greedy (default) or dynamic
#   cidder protein clustering (CD-HIT) of the proteins predicted by Pyrodigal (-a): representatives are added until
#          they sample the pangenome of the taxon (CiDDER defaults). If CiDDER cannot read the proteomes (CiDDER
#          1.3.6 fails without genome files), the genomes are given instead and CiDDER predicts their genes itself
# --filter-mge makes skDER/CiDDER mask the mobile genetic elements (prophages and plasmids, predicted with PhiSpy)
# before comparing the genomes.
#
# All the results go to one folder (-o, OUTPUT/GENOME-DEREPLICATION), not to the leaves of the hierarchy:
#   dereplication-summary.tsv     one row per lineage: TAXONOMIC_LINEAGE (Domain;...;Species), TAXONOMIC_LEAF_RANK,
#                                 GENOMES_BEFORE_DEREPLICATION, GENOMES_AFTER_DEREPLICATION, REPRESENTATIVE_GENOMES
#                                 (accessions, ";"), DEREPLICATION_SOFTWARE (skder | cidder | NONE for a single
#                                 genome), SKDER_MODE, SKDER_ANI_CUTOFF, CIDDER_INPUT (proteomes | genomes),
#                                 MGE_FILTERED (yes | no), MGE_SOFTWARE (PhiSpy), DEREPLICATION_STATUS
#                                 (DEREPLICATED | SINGLE_GENOME | FAILED)
#   representative-genomes.tsv    ASSEMBLY_ACCESSION  TAXONOMIC_LINEAGE of every representative genome
#   genome-dereplication.tsv      (or -t) one row per genome: ASSEMBLY_ACCESSION  REPRESENTATIVE_GENOME
#                                 (yes | no | NOT_DEREPLICATED)  LEAF_GENOME_COUNT  LEAF_REPRESENTATIVE_COUNT
#   RUNS/<Domain>/.../<Species>/  raw output of skDER/CiDDER of each lineage (also the cache of the re-execution)
# Each leaf of the hierarchy also gets GENOME-METADATA/genome-metadata.tsv (its genomes + REPRESENTATIVE_GENOME).
#
# Re-execution: the result of a lineage is reused when its list of genomes and the parameters did not change.
#
# Usage: dereplicate-genomes.sh -d <hierarchy_dir> -m <genome_metadata.tsv> -o <dereplication_dir> [-t <table>]
#            [-p <threads>] [--memory <GB>] [--method skder|cidder] [-a <genes_dir>] [--mode greedy|dynamic] [--ani <pct>]
#            [--filter-mge]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -d <hierarchy_dir> -m <genome_metadata.tsv> -o <dereplication_dir> [-t <table>] [options]

  -d <dir>            Root of the taxonomic hierarchy (with leaf-directories.tsv)
  -m <file>           Genome metadata table (copied, per leaf, to <leaf>/GENOME-METADATA)
  -o <dir>            Folder of every dereplication result (tables and RUNS/)
  -t <file>           Per-genome result table (default: <dereplication_dir>/genome-dereplication.tsv)
  -p <n>              CPUs (default: 4)
  --memory <GB>       Maximum RAM of skDER/CiDDER (default: no limit)
  --method <method>   skder (default; ANI) or cidder (protein clustering of the predicted proteins)
  -a <dir>            Folder with the predicted proteins, <ACCESSION>.faa (used by cidder)
  --mode <mode>       skDER dereplication mode: greedy (default) or dynamic
  --ani <pct>         skDER ANI (% identity) cutoff: genomes above it are redundant (default: 99.5)
  --filter-mge        Mask mobile genetic elements (PhiSpy) before comparing the genomes
EOF
    exit "${1:-1}"
}

HIERARCHY_DIR=""
METADATA=""
TABLE=""
DEREP_OUT=""
THREADS=4
MEMORY_GB=0   # --max-memory of skDER/CiDDER (0 = no limit)
GENES_DIR=""
DEREP_METHOD="skder"
SKDER_MODE="greedy"
SKDER_PERCENT_IDENTITY="99.5"
SKDER_ALIGNED_FRACTION="90.0"
SKDER_EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d) HIERARCHY_DIR="${2%/}"; shift 2 ;;
        -m) METADATA="$2"; shift 2 ;;
        -t) TABLE="$2"; shift 2 ;;
        -o) DEREP_OUT="${2%/}"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        --memory) MEMORY_GB="$2"; shift 2 ;;
        -a) GENES_DIR="${2%/}"; shift 2 ;;
        --method) DEREP_METHOD="${2,,}"; shift 2 ;;
        --mode) SKDER_MODE="${2,,}"; shift 2 ;;
        --ani) SKDER_PERCENT_IDENTITY="$2"; shift 2 ;;
        --filter-mge) SKDER_EXTRA_ARGS="--filter-mge"; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${HIERARCHY_DIR}" || -z "${METADATA}" || -z "${DEREP_OUT}" ]]; then
    log "error" "Options -d, -m and -o are required."
    usage
fi
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "The number of threads must be a positive integer (got '${THREADS}')."
[[ "${MEMORY_GB}" =~ ^[0-9]+$ ]] || die "--memory must be a number of GB (got '${MEMORY_GB}')."
[[ "${SKDER_MODE}" == "greedy" || "${SKDER_MODE}" == "dynamic" ]] || die "Invalid skDER mode '${SKDER_MODE}'. Must be greedy or dynamic."
[[ "${DEREP_METHOD}" == "skder" || "${DEREP_METHOD}" == "cidder" ]] || die "Invalid dereplication method '${DEREP_METHOD}'. Must be skder or cidder."
if [[ "${DEREP_METHOD}" == "cidder" && ( -z "${GENES_DIR}" || ! -d "${GENES_DIR}" ) ]]; then
    die "--method cidder needs the folder of the predicted proteins (-a)."
fi
if [[ ! "${SKDER_PERCENT_IDENTITY}" =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v v="${SKDER_PERCENT_IDENTITY}" 'BEGIN { exit !(v + 0 > 0 && v + 0 <= 100) }'; then
    die "Invalid ANI cutoff '${SKDER_PERCENT_IDENTITY}'. It must be a percentage greater than 0 and at most 100."
fi

LEAF_TABLE="${HIERARCHY_DIR}/leaf-directories.tsv"
[[ -f "${LEAF_TABLE}" ]] || die "Leaf table not found: ${LEAF_TABLE}. Run build-taxonomic-dir-hierarchy.sh first."
[[ -f "${METADATA}" ]] || die "Metadata table not found: ${METADATA}"

mkdir -p "${DEREP_OUT}"
setup_module_log "${DEREP_OUT}/dereplicate-genomes.log"

LIST_NAME="genomes-files-symlink.txt"
RUNS_DIR="${DEREP_OUT}/RUNS"
SUMMARY_TABLE="${DEREP_OUT}/dereplication-summary.tsv"
REPS_TABLE="${DEREP_OUT}/representative-genomes.tsv"
[[ -n "${TABLE}" ]] || TABLE="${DEREP_OUT}/genome-dereplication.tsv"
HIERARCHY_ABS="$(absolute_path "${HIERARCHY_DIR}")"
MGE_SOFTWARE="NULL"
[[ -z "${SKDER_EXTRA_ARGS}" ]] || MGE_SOFTWARE="PhiSpy"

make_temp_dir WORK_DIR

# Path of a leaf below the hierarchy (Bacteria/.../Escherichia_coli); the lineage is the same with ";"
leaf_relative_path() {
    local REL="${1#"${HIERARCHY_ABS}"/}"
    echo "${REL#"${HIERARCHY_DIR}"/}"
}

# Runs CiDDER on the predicted proteins of the genomes of a leaf; if CiDDER cannot use them (CiDDER 1.3.6 stops when
# no genome file is given), runs it again on the genome files (CiDDER then predicts the genes itself). The input
# used is written to <run_dir>/.cidder-input (proteomes or genomes).
# Usage: run_cidder <threads> <leaf_dir> <run_dir> <genome files...>
run_cidder() {
    local THREADS="$1" LEAF="$2" RUN_DIR="$3"
    local -a PROTEOMES=()
    local ACC
    shift 3
    while IFS=$'\t' read -r ACC _; do
        PROTEOMES+=("${GENES_DIR}/${ACC}.faa")
    done < "${LEAF}/${LIST_NAME}"
    # shellcheck disable=SC2086
    if cidder --proteomes "${PROTEOMES[@]}" --output-directory "${RUN_DIR}/CIDDER-RESULTS" --threads "${THREADS}" \
            --max-memory "${MEMORY_GB}" --symlink ${SKDER_EXTRA_ARGS} < /dev/null > "${RUN_DIR}/cidder-console.log" 2>&1; then
        echo "proteomes" > "${RUN_DIR}/.cidder-input"
        return 0
    fi
    log "warning" "CiDDER could not use the predicted proteins of $(leaf_relative_path "${LEAF}") (see ${RUN_DIR}/cidder-console.log); running it on the genome files (CiDDER predicts their genes)."
    rm -rf "${RUN_DIR}/CIDDER-RESULTS"
    echo "genomes" > "${RUN_DIR}/.cidder-input"
    # shellcheck disable=SC2086
    cidder --genomes "$@" --output-directory "${RUN_DIR}/CIDDER-RESULTS" --threads "${THREADS}" \
        --max-memory "${MEMORY_GB}" --symlink ${SKDER_EXTRA_ARGS} < /dev/null >> "${RUN_DIR}/cidder-console.log" 2>&1
}

# Dereplicates one leaf with skDER or CiDDER (one leaf at a time, with every CPU), in RUNS/<leaf path>/. Writes
# representative-genomes.txt there; on failure the folder gets a FAILED marker instead.
# Usage: dereplicate_leaf <threads> <leaf_dir>
dereplicate_leaf() {
    local THREADS="$1" LEAF="$2"
    local RUN_DIR
    local PARAMS="method=${DEREP_METHOD} extra=${SKDER_EXTRA_ARGS}"
    local INPUT_LIST REPS_DIR RESULTS_DIR STATUS
    local -a GENOMES

    RUN_DIR="${RUNS_DIR}/$(leaf_relative_path "${LEAF}")"
    if [[ "${DEREP_METHOD}" == "skder" ]]; then
        PARAMS="${PARAMS} mode=${SKDER_MODE} identity=${SKDER_PERCENT_IDENTITY} aligned_fraction=${SKDER_ALIGNED_FRACTION}"
        RESULTS_DIR="${RUN_DIR}/SKDER-RESULTS"
    else
        RESULTS_DIR="${RUN_DIR}/CIDDER-RESULTS"
    fi
    INPUT_LIST="$(cut -f2 "${LEAF}/${LIST_NAME}" | LC_ALL=C sort)"

    # Same genomes and parameters as the previous run: reuse
    if [[ -s "${RUN_DIR}/representative-genomes.txt" && -f "${RUN_DIR}/.input-genomes.txt" && -f "${RUN_DIR}/.parameters.txt" ]] \
       && [[ "$(cat "${RUN_DIR}/.input-genomes.txt")" == "${INPUT_LIST}" && "$(cat "${RUN_DIR}/.parameters.txt")" == "${PARAMS}" ]]; then
        return 0
    fi

    rm -rf "${RUN_DIR}"
    mkdir -p "${RUN_DIR}"
    mapfile -t GENOMES <<< "${INPUT_LIST}"

    STATUS=0
    if [[ "${DEREP_METHOD}" == "skder" ]]; then
        # shellcheck disable=SC2086
        skder --genomes "${GENOMES[@]}" \
              --output-directory "${RESULTS_DIR}" \
              --dereplication-mode "${SKDER_MODE}" \
              --percent-identity-cutoff "${SKDER_PERCENT_IDENTITY}" \
              --aligned-fraction-cutoff "${SKDER_ALIGNED_FRACTION}" \
              --threads "${THREADS}" \
              --max-memory "${MEMORY_GB}" \
              --symlink \
              ${SKDER_EXTRA_ARGS} < /dev/null > "${RUN_DIR}/skder-console.log" 2>&1 || STATUS=$?
    else
        run_cidder "${THREADS}" "${LEAF}" "${RUN_DIR}" "${GENOMES[@]}" || STATUS=$?
    fi
    if [[ "${STATUS}" -ne 0 ]]; then
        touch "${RUN_DIR}/FAILED"
        log "warning" "${DEREP_METHOD} failed in $(leaf_relative_path "${LEAF}") (see ${RUN_DIR}/${DEREP_METHOD}-console.log); its genomes are marked NOT_DEREPLICATED."
        return 0
    fi

    # Representatives = files linked in the output folder (their names start with the accession)
    REPS_DIR="$(find "${RESULTS_DIR}" -maxdepth 1 -type d -name 'Dereplicated_Representative_Genomes' | head -n 1)"
    if [[ -n "${REPS_DIR}" ]]; then
        find "${REPS_DIR}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -printf '%f\n' \
            | grep -oE '^GC[AF]_[0-9]+[.][0-9]+' | LC_ALL=C sort -u > "${RUN_DIR}/representative-genomes.txt" || true
    fi
    if [[ ! -s "${RUN_DIR}/representative-genomes.txt" ]]; then
        touch "${RUN_DIR}/FAILED"
        log "warning" "${DEREP_METHOD} wrote no representative genome for $(leaf_relative_path "${LEAF}"); its genomes are marked NOT_DEREPLICATED."
        return 0
    fi
    printf '%s\n' "${INPUT_LIST}" > "${RUN_DIR}/.input-genomes.txt"
    printf '%s\n' "${PARAMS}" > "${RUN_DIR}/.parameters.txt"
}

# 0. Clean-up: dereplication folders inside the hierarchy (written by older versions of Prokfinder) and runs of
#    lineages that no longer exist
find "${HIERARCHY_DIR}" -mindepth 2 -type d -name "GENOME-DEREPLICATION" -prune -exec rm -rf {} +
if [[ -d "${RUNS_DIR}" ]]; then
    tail -n +2 "${LEAF_TABLE}" | cut -f1 | while IFS= read -r LEAF; do
        echo "${RUNS_DIR}/$(leaf_relative_path "${LEAF}")"
    done | LC_ALL=C sort > "${WORK_DIR}/current_runs.txt"
    find "${RUNS_DIR}" -type f -name ".parameters.txt" -printf '%h\n' | LC_ALL=C sort \
        | LC_ALL=C comm -23 - "${WORK_DIR}/current_runs.txt" | while IFS= read -r OLD; do rm -rf -- "${OLD}"; done
    find "${RUNS_DIR}" -mindepth 1 -type d -empty -delete
fi

# 1. Leaves with a single genome need no dereplication; the others are dereplicated one at a time with every CPU
tail -n +2 "${LEAF_TABLE}" | awk -F'\t' '$3 > 1 { print $1 }' > "${WORK_DIR}/multi_leaves.txt"
N_MULTI="$(count_lines "${WORK_DIR}/multi_leaves.txt")"
log "message" "Dereplication: $(( $(count_records "${LEAF_TABLE}") - N_MULTI )) lineage(s) with a single genome; ${N_MULTI} with several genomes (${DEREP_METHOD})."

if [[ "${N_MULTI}" -gt 0 ]]; then
    require_programs "${DEREP_METHOD}"
    if [[ "${DEREP_METHOD}" == "skder" ]]; then
        log "message" "skDER: $(skder --version 2>&1 | head -n 1); mode ${SKDER_MODE}, identity ${SKDER_PERCENT_IDENTITY}%, aligned fraction ${SKDER_ALIGNED_FRACTION}%; ${THREADS} CPU(s)."
    else
        log "message" "CiDDER: $(cidder --version 2>&1 | head -n 1); predicted proteins from ${GENES_DIR}; CiDDER default parameters; ${THREADS} CPU(s)."
    fi
    while IFS= read -r LEAF; do
        dereplicate_leaf "${THREADS}" "${LEAF}"
    done < "${WORK_DIR}/multi_leaves.txt"
fi

# 2. Tables: per lineage (summary), per genome, and the list of representatives
#    <leaf> <lineage> <rank> <genomes> <status> <representatives (;)> -> written to the summary and used by the others
{
    printf 'TAXONOMIC_LINEAGE\tTAXONOMIC_LEAF_RANK\tGENOMES_BEFORE_DEREPLICATION\tGENOMES_AFTER_DEREPLICATION\tREPRESENTATIVE_GENOMES\t'
    printf 'DEREPLICATION_SOFTWARE\tSKDER_MODE\tSKDER_ANI_CUTOFF\tCIDDER_INPUT\tMGE_FILTERED\tMGE_SOFTWARE\tDEREPLICATION_STATUS\n'
    tail -n +2 "${LEAF_TABLE}" | while IFS=$'\t' read -r LEAF RANK COUNT; do
        REL="$(leaf_relative_path "${LEAF}")"
        LINEAGE="${REL//\//;}"
        RUN_DIR="${RUNS_DIR}/${REL}"
        if [[ "${COUNT}" -eq 1 ]]; then
            printf '%s\t%s\t1\t1\t%s\tNONE\tNULL\tNULL\tNULL\tno\tNULL\tSINGLE_GENOME\n' "${LINEAGE}" "${RANK}" "$(cut -f1 "${LEAF}/${LIST_NAME}")"
            continue
        fi
        SKDER_COLS="NULL\tNULL"
        [[ "${DEREP_METHOD}" != "skder" ]] || SKDER_COLS="${SKDER_MODE}\t${SKDER_PERCENT_IDENTITY}"
        CIDDER_INPUT="NULL"
        [[ "${DEREP_METHOD}" != "cidder" ]] || CIDDER_INPUT="$(cat "${RUN_DIR}/.cidder-input" 2>/dev/null || echo NULL)"
        MGE="no"
        [[ -z "${SKDER_EXTRA_ARGS}" ]] || MGE="yes"
        if [[ -f "${RUN_DIR}/FAILED" || ! -s "${RUN_DIR}/representative-genomes.txt" ]]; then
            printf "%s\t%s\t%s\tNULL\tNULL\t%s\t${SKDER_COLS}\t%s\t%s\t%s\tFAILED\n" \
                "${LINEAGE}" "${RANK}" "${COUNT}" "${DEREP_METHOD}" "${CIDDER_INPUT}" "${MGE}" "${MGE_SOFTWARE}"
        else
            printf "%s\t%s\t%s\t%s\t%s\t%s\t${SKDER_COLS}\t%s\t%s\t%s\tDEREPLICATED\n" \
                "${LINEAGE}" "${RANK}" "${COUNT}" "$(count_lines "${RUN_DIR}/representative-genomes.txt")" \
                "$(paste -sd ';' "${RUN_DIR}/representative-genomes.txt")" "${DEREP_METHOD}" "${CIDDER_INPUT}" \
                "${MGE}" "${MGE_SOFTWARE}"
        fi
    done
} > "${SUMMARY_TABLE}"

# Per genome (also added to genome-metadata.tsv) and list of representatives
{
    printf 'ASSEMBLY_ACCESSION\tREPRESENTATIVE_GENOME\tLEAF_GENOME_COUNT\tLEAF_REPRESENTATIVE_COUNT\n'
    tail -n +2 "${LEAF_TABLE}" | while IFS=$'\t' read -r LEAF _ COUNT; do
        LINEAGE="$(leaf_relative_path "${LEAF}")"
        LINEAGE="${LINEAGE//\//;}"
        awk -F'\t' -v OFS='\t' -v lineage="${LINEAGE}" -v n="${COUNT}" '
            NR == FNR { if ($1 == lineage) { st = $12; nr = $4; m = split($5, r, ";"); for (i = 1; i <= m; i++) rep[r[i]] = 1 } next }
            {
                if (st == "FAILED") print $1, "NOT_DEREPLICATED", n, "NULL"
                else print $1, (($1 in rep) ? "yes" : "no"), n, nr
            }' "${SUMMARY_TABLE}" "${LEAF}/${LIST_NAME}"
    done | LC_ALL=C sort -t $'\t' -k1,1
} > "${TABLE}"

{
    printf 'ASSEMBLY_ACCESSION\tTAXONOMIC_LINEAGE\n'
    awk -F'\t' -v OFS='\t' 'NR > 1 && $12 != "FAILED" { n = split($5, r, ";"); for (i = 1; i <= n; i++) print r[i], $1 }' \
        "${SUMMARY_TABLE}" | LC_ALL=C sort -t $'\t' -k1,1
} > "${REPS_TABLE}"

# 3. Metadata of each leaf of the hierarchy (rows of its genomes + REPRESENTATIVE_GENOME)
tail -n +2 "${LEAF_TABLE}" | cut -f1 | while IFS= read -r LEAF; do
    mkdir -p "${LEAF}/GENOME-METADATA"
done

# (add_metadata_columns replaces the dereplication columns if the metadata already has them from a previous run)
cp "${METADATA}" "${WORK_DIR}/metadata.tsv"
add_metadata_columns "${WORK_DIR}/metadata.tsv" "${TABLE}"

awk -F'\t' -v OFS='\t' -v placement="${HIERARCHY_DIR}/genome-dir-hierarchy.tsv" '
    BEGIN {
        while ((getline line < placement) > 0) { split(line, f, "\t"); if (f[1] != "ASSEMBLY_ACCESSION") leaf[f[1]] = f[3] }
        close(placement)
    }
    NR == 1 { next }
    ($1 in leaf) { print leaf[$1], $0 }
' "${WORK_DIR}/metadata.tsv" | LC_ALL=C sort -t $'\t' -k1,1 | awk -F'\t' -v header="$(head -n 1 "${WORK_DIR}/metadata.tsv")" '
    $1 != current {
        if (current != "") close(out)
        current = $1
        out = current "/GENOME-METADATA/genome-metadata.tsv"
        print header > out
    }
    { sub(/^[^\t]*\t/, ""); print > out }
'

N_REP="$(awk -F'\t' 'NR > 1 && $2 == "yes"' "${TABLE}" | wc -l)"
N_FAILED="$(awk -F'\t' 'NR > 1 && $2 == "NOT_DEREPLICATED"' "${TABLE}" | wc -l)"
log "message" "Dereplication complete: ${N_REP} representative genome(s) out of $(count_records "${TABLE}"), in $(count_records "${SUMMARY_TABLE}") lineage(s)."
log "message" "  Per lineage: ${SUMMARY_TABLE}"
log "message" "  Representatives: ${REPS_TABLE}"
if [[ "${N_FAILED}" -gt 0 ]]; then
    log "warning" "${N_FAILED} genome(s) could not be dereplicated (NOT_DEREPLICATED); run Prokfinder again to retry."
fi
