#!/usr/bin/env bash
#
# dereplicate-genomes-skder.sh - Prokfinder module 7: dereplication of the genomes of each taxon
#
# skDER is run inside each leaf directory of the taxonomic hierarchy, i.e. only genomes with the same lineage
# (same species, or the lowest known rank) are compared. Leaves with a single genome need no dereplication: the
# genome is its own representative. Depending on the dereplication parameters, a leaf may keep more than one
# representative genome when there is divergence/variability inside the taxon.
#
# Re-execution: the result of a leaf is reused when its list of genomes and the skDER parameters did not change.
#
# Usage: dereplicate-genomes-skder.sh -d <hierarchy_dir> -m <genome_metadata.tsv> -t <table> [-p <threads>]
#
# Output:
#   <leaf>/GENOME-DEREPLICATION/              skDER output + representative-genomes.txt
#   <leaf>/GENOME-METADATA/genome-metadata.tsv   metadata of the genomes of the leaf + REPRESENTATIVE_GENOME
#   <table>: ASSEMBLY_ACCESSION  REPRESENTATIVE_GENOME (yes | no | NOT_DEREPLICATED)  LEAF_GENOME_COUNT  LEAF_REPRESENTATIVE_COUNT

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -d <hierarchy_dir> -m <genome_metadata.tsv> -t <table> [-p <threads>]

  -d <dir>    Root of the taxonomic hierarchy (with leaf-directories.tsv)
  -m <file>   Genome metadata table (copied, per leaf, to <leaf>/GENOME-METADATA)
  -t <file>   Result table (one row per genome)
  -p <n>      Total CPUs (default: number of CPUs)

Environment variables:
  SKDER_MODE (greedy), SKDER_PERCENT_IDENTITY (99.5, the skDER default), SKDER_ALIGNED_FRACTION (90.0),
  SKDER_THREADS_PER_JOB (4: CPUs of each skDER run; several leaves run at the same time),
  SKDER_EXTRA_ARGS (extra options passed to skder)
EOF
    exit "${1:-1}"
}

HIERARCHY_DIR=""
METADATA=""
TABLE=""
THREADS="$(nproc 2>/dev/null || echo 4)"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d) HIERARCHY_DIR="${2%/}"; shift 2 ;;
        -m) METADATA="$2"; shift 2 ;;
        -t) TABLE="$2"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${HIERARCHY_DIR}" || -z "${METADATA}" || -z "${TABLE}" ]]; then
    log "error" "Options -d, -m and -t are required."
    usage
fi
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "The number of threads must be a positive integer (got '${THREADS}')."
[[ "${SKDER_THREADS_PER_JOB:-4}" =~ ^[1-9][0-9]*$ ]] || die "SKDER_THREADS_PER_JOB must be a positive integer."

LEAF_TABLE="${HIERARCHY_DIR}/leaf-directories.tsv"
[[ -f "${LEAF_TABLE}" ]] || die "Leaf table not found: ${LEAF_TABLE}. Run build-taxonomic-dir-hierarchy.sh first."
[[ -f "${METADATA}" ]] || die "Metadata table not found: ${METADATA}"

setup_module_log "${HIERARCHY_DIR}/dereplicate-genomes-skder.log"

export SKDER_MODE="${SKDER_MODE:-greedy}"
export SKDER_PERCENT_IDENTITY="${SKDER_PERCENT_IDENTITY:-99.5}"
export SKDER_ALIGNED_FRACTION="${SKDER_ALIGNED_FRACTION:-90.0}"
export SKDER_EXTRA_ARGS="${SKDER_EXTRA_ARGS:-}"
SKDER_THREADS_PER_JOB="${SKDER_THREADS_PER_JOB:-4}"
LIST_NAME="genomes-files-symlink.txt"

make_temp_dir WORK_DIR

# Dereplicates one leaf with skDER (runs in parallel through xargs). Writes
# <leaf>/GENOME-DEREPLICATION/representative-genomes.txt; on failure the folder gets a FAILED marker instead.
# Usage: dereplicate_leaf <threads> <leaf_dir>
dereplicate_leaf() {
    local THREADS="$1" LEAF="$2"
    local DEREP_DIR="${LEAF}/GENOME-DEREPLICATION"
    local PARAMS="mode=${SKDER_MODE} identity=${SKDER_PERCENT_IDENTITY} aligned_fraction=${SKDER_ALIGNED_FRACTION} extra=${SKDER_EXTRA_ARGS}"
    local INPUT_LIST REPS_DIR
    local -a GENOMES

    INPUT_LIST="$(cut -f2 "${LEAF}/genomes-files-symlink.txt" | LC_ALL=C sort)"

    # Same genomes and parameters as the previous run: reuse
    if [[ -s "${DEREP_DIR}/representative-genomes.txt" && -f "${DEREP_DIR}/.input-genomes.txt" && -f "${DEREP_DIR}/.parameters.txt" ]] \
       && [[ "$(cat "${DEREP_DIR}/.input-genomes.txt")" == "${INPUT_LIST}" && "$(cat "${DEREP_DIR}/.parameters.txt")" == "${PARAMS}" ]]; then
        return 0
    fi

    rm -rf "${DEREP_DIR}"
    mkdir -p "${DEREP_DIR}"
    mapfile -t GENOMES <<< "${INPUT_LIST}"

    # shellcheck disable=SC2086
    if ! skder --genomes "${GENOMES[@]}" \
               --output-directory "${DEREP_DIR}/SKDER-RESULTS" \
               --dereplication-mode "${SKDER_MODE}" \
               --percent-identity-cutoff "${SKDER_PERCENT_IDENTITY}" \
               --aligned-fraction-cutoff "${SKDER_ALIGNED_FRACTION}" \
               --threads "${THREADS}" \
               --symlink \
               ${SKDER_EXTRA_ARGS} > "${DEREP_DIR}/skder-console.log" 2>&1; then
        touch "${DEREP_DIR}/FAILED"
        log "warning" "skDER failed in ${LEAF} (see ${DEREP_DIR}/skder-console.log); its genomes are marked NOT_DEREPLICATED."
        return 0
    fi

    # Representatives = files linked in skDER's output folder (their names start with the accession)
    REPS_DIR="$(find "${DEREP_DIR}/SKDER-RESULTS" -maxdepth 1 -type d -name 'Dereplicated_Representative_Genomes' | head -n 1)"
    if [[ -z "${REPS_DIR}" ]]; then
        touch "${DEREP_DIR}/FAILED"
        log "warning" "skDER wrote no representative genome in ${LEAF}; its genomes are marked NOT_DEREPLICATED."
        return 0
    fi
    find "${REPS_DIR}" -mindepth 1 -maxdepth 1 \( -type f -o -type l \) -printf '%f\n' \
        | grep -oE '^GC[AF]_[0-9]+[.][0-9]+' | LC_ALL=C sort -u > "${DEREP_DIR}/representative-genomes.txt" || true

    if [[ ! -s "${DEREP_DIR}/representative-genomes.txt" ]]; then
        touch "${DEREP_DIR}/FAILED"
        log "warning" "No representative genome could be read from ${REPS_DIR}; the genomes of ${LEAF} are marked NOT_DEREPLICATED."
        return 0
    fi
    printf '%s\n' "${INPUT_LIST}" > "${DEREP_DIR}/.input-genomes.txt"
    printf '%s\n' "${PARAMS}" > "${DEREP_DIR}/.parameters.txt"
}

# 1. Leaves with a single genome: the genome is the representative (no skDER run)
tail -n +2 "${LEAF_TABLE}" | awk -F'\t' '$3 == 1 { print $1 }' > "${WORK_DIR}/single_leaves.txt"
tail -n +2 "${LEAF_TABLE}" | awk -F'\t' '$3 > 1 { print $1 }' > "${WORK_DIR}/multi_leaves.txt"

while IFS= read -r LEAF; do
    DEREP_DIR="${LEAF}/GENOME-DEREPLICATION"
    rm -rf "${DEREP_DIR}"
    mkdir -p "${DEREP_DIR}"
    cut -f1 "${LEAF}/${LIST_NAME}" > "${DEREP_DIR}/representative-genomes.txt"
done < "${WORK_DIR}/single_leaves.txt"

N_MULTI="$(count_lines "${WORK_DIR}/multi_leaves.txt")"
log "message" "Dereplication: $(count_lines "${WORK_DIR}/single_leaves.txt") leaf directories with a single genome; ${N_MULTI} with several genomes (skDER)."

# 2. Leaves with several genomes: skDER, several leaves at the same time
if [[ "${N_MULTI}" -gt 0 ]]; then
    require_programs skder
    JOB_THREADS=$(( SKDER_THREADS_PER_JOB < THREADS ? SKDER_THREADS_PER_JOB : THREADS ))
    JOBS=$(( THREADS / JOB_THREADS ))
    log "message" "skDER: $(skder --version 2>&1 | head -n 1); mode ${SKDER_MODE}, identity ${SKDER_PERCENT_IDENTITY}%, aligned fraction ${SKDER_ALIGNED_FRACTION}%; ${JOBS} leaf(s) at a time with ${JOB_THREADS} CPU(s) each."
    export -f log dereplicate_leaf
    xargs -a "${WORK_DIR}/multi_leaves.txt" -d '\n' -n 1 -P "${JOBS}" bash -c 'dereplicate_leaf "$@"' _ "${JOB_THREADS}" \
        || die "Dereplication stopped unexpectedly."
fi

# 3. Result table
{
    printf 'ASSEMBLY_ACCESSION\tREPRESENTATIVE_GENOME\tLEAF_GENOME_COUNT\tLEAF_REPRESENTATIVE_COUNT\n'
    tail -n +2 "${LEAF_TABLE}" | while IFS=$'\t' read -r LEAF _ COUNT; do
        REPS_FILE="${LEAF}/GENOME-DEREPLICATION/representative-genomes.txt"
        if [[ -f "${LEAF}/GENOME-DEREPLICATION/FAILED" || ! -s "${REPS_FILE}" ]]; then
            awk -F'\t' -v OFS='\t' -v n="${COUNT}" '{ print $1, "NOT_DEREPLICATED", n, "NULL" }' "${LEAF}/${LIST_NAME}"
        else
            awk -F'\t' -v OFS='\t' -v n="${COUNT}" 'NR == FNR { rep[$1] = 1; nrep++; next }
                { print $1, (($1 in rep) ? "yes" : "no"), n, nrep }' "${REPS_FILE}" "${LEAF}/${LIST_NAME}"
        fi
    done | LC_ALL=C sort -t $'\t' -k1,1
} > "${TABLE}"

# 4. Metadata of each leaf (rows of its genomes + REPRESENTATIVE_GENOME)
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
log "message" "Dereplication complete: ${N_REP} representative genome(s) out of $(count_records "${TABLE}"). Table: ${TABLE}"
if [[ "${N_FAILED}" -gt 0 ]]; then
    log "warning" "${N_FAILED} genome(s) could not be dereplicated (NOT_DEREPLICATED); run Prokfinder again to retry."
fi
