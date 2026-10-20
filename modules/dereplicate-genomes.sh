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
# Re-execution: the result of a leaf is reused when its list of genomes and the parameters did not change.
#
# Usage: dereplicate-genomes.sh -d <hierarchy_dir> -m <genome_metadata.tsv> -t <table> [-p <threads>]
#            [--method skder|cidder] [-a <genes_dir>] [--mode greedy|dynamic] [--ani <pct>] [--filter-mge]
#
# Output:
#   <leaf>/GENOME-DEREPLICATION/              skDER/CiDDER output + representative-genomes.txt
#   <leaf>/GENOME-METADATA/genome-metadata.tsv   metadata of the genomes of the leaf + REPRESENTATIVE_GENOME
#   <table>: ASSEMBLY_ACCESSION  REPRESENTATIVE_GENOME (yes | no | NOT_DEREPLICATED)  LEAF_GENOME_COUNT  LEAF_REPRESENTATIVE_COUNT

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -d <hierarchy_dir> -m <genome_metadata.tsv> -t <table> [-p <threads>] [--mode <mode>] [--filter-mge]

  -d <dir>        Root of the taxonomic hierarchy (with leaf-directories.tsv)
  -m <file>       Genome metadata table (copied, per leaf, to <leaf>/GENOME-METADATA)
  -t <file>       Result table (one row per genome)
  -p <n>              CPUs (default: 4)
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
THREADS=4
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
        -p) THREADS="$2"; shift 2 ;;
        -a) GENES_DIR="${2%/}"; shift 2 ;;
        --method) DEREP_METHOD="${2,,}"; shift 2 ;;
        --mode) SKDER_MODE="${2,,}"; shift 2 ;;
        --ani) SKDER_PERCENT_IDENTITY="$2"; shift 2 ;;
        --filter-mge) SKDER_EXTRA_ARGS="--filter-mge"; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${HIERARCHY_DIR}" || -z "${METADATA}" || -z "${TABLE}" ]]; then
    log "error" "Options -d, -m and -t are required."
    usage
fi
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "The number of threads must be a positive integer (got '${THREADS}')."
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

setup_module_log "${HIERARCHY_DIR}/dereplicate-genomes.log"

LIST_NAME="genomes-files-symlink.txt"

make_temp_dir WORK_DIR

# Runs CiDDER on the predicted proteins of the genomes of a leaf; if CiDDER cannot use them (CiDDER 1.3.6 stops when
# no genome file is given), runs it again on the genome files (CiDDER then predicts the genes itself)
# Usage: run_cidder <threads> <leaf_dir> <derep_dir> <genome files...>
run_cidder() {
    local THREADS="$1" LEAF="$2" DEREP_DIR="$3"
    local -a PROTEOMES=()
    local ACC
    shift 3
    while IFS=$'\t' read -r ACC _; do
        PROTEOMES+=("${GENES_DIR}/${ACC}.faa")
    done < "${LEAF}/${LIST_NAME}"
    # shellcheck disable=SC2086
    if cidder --proteomes "${PROTEOMES[@]}" --output-directory "${DEREP_DIR}/CIDDER-RESULTS" --threads "${THREADS}" \
            --symlink ${SKDER_EXTRA_ARGS} < /dev/null > "${DEREP_DIR}/cidder-console.log" 2>&1; then
        return 0
    fi
    log "warning" "CiDDER could not use the predicted proteins in ${LEAF} (see ${DEREP_DIR}/cidder-console.log); running it on the genome files (CiDDER predicts their genes)."
    rm -rf "${DEREP_DIR}/CIDDER-RESULTS"
    # shellcheck disable=SC2086
    cidder --genomes "$@" --output-directory "${DEREP_DIR}/CIDDER-RESULTS" --threads "${THREADS}" \
        --symlink ${SKDER_EXTRA_ARGS} < /dev/null >> "${DEREP_DIR}/cidder-console.log" 2>&1
}

# Dereplicates one leaf with skDER or CiDDER (one leaf at a time, with every CPU). Writes
# <leaf>/GENOME-DEREPLICATION/representative-genomes.txt; on failure the folder gets a FAILED marker instead.
# Usage: dereplicate_leaf <threads> <leaf_dir>
dereplicate_leaf() {
    local THREADS="$1" LEAF="$2"
    local DEREP_DIR="${LEAF}/GENOME-DEREPLICATION"
    local PARAMS="method=${DEREP_METHOD} extra=${SKDER_EXTRA_ARGS}"
    local INPUT_LIST REPS_DIR RESULTS_DIR STATUS
    local -a GENOMES

    if [[ "${DEREP_METHOD}" == "skder" ]]; then
        PARAMS="${PARAMS} mode=${SKDER_MODE} identity=${SKDER_PERCENT_IDENTITY} aligned_fraction=${SKDER_ALIGNED_FRACTION}"
        RESULTS_DIR="${DEREP_DIR}/SKDER-RESULTS"
    else
        RESULTS_DIR="${DEREP_DIR}/CIDDER-RESULTS"
    fi
    INPUT_LIST="$(cut -f2 "${LEAF}/genomes-files-symlink.txt" | LC_ALL=C sort)"

    # Same genomes and parameters as the previous run: reuse
    if [[ -s "${DEREP_DIR}/representative-genomes.txt" && -f "${DEREP_DIR}/.input-genomes.txt" && -f "${DEREP_DIR}/.parameters.txt" ]] \
       && [[ "$(cat "${DEREP_DIR}/.input-genomes.txt")" == "${INPUT_LIST}" && "$(cat "${DEREP_DIR}/.parameters.txt")" == "${PARAMS}" ]]; then
        return 0
    fi

    rm -rf "${DEREP_DIR}"
    mkdir -p "${DEREP_DIR}"
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
              --symlink \
              ${SKDER_EXTRA_ARGS} < /dev/null > "${DEREP_DIR}/skder-console.log" 2>&1 || STATUS=$?
    else
        run_cidder "${THREADS}" "${LEAF}" "${DEREP_DIR}" "${GENOMES[@]}" || STATUS=$?
    fi
    if [[ "${STATUS}" -ne 0 ]]; then
        touch "${DEREP_DIR}/FAILED"
        log "warning" "${DEREP_METHOD} failed in ${LEAF} (see ${DEREP_DIR}/${DEREP_METHOD}-console.log); its genomes are marked NOT_DEREPLICATED."
        return 0
    fi

    # Representatives = files linked in the output folder (their names start with the accession)
    REPS_DIR="$(find "${RESULTS_DIR}" -maxdepth 1 -type d -name 'Dereplicated_Representative_Genomes' | head -n 1)"
    if [[ -z "${REPS_DIR}" ]]; then
        touch "${DEREP_DIR}/FAILED"
        log "warning" "${DEREP_METHOD} wrote no representative genome in ${LEAF}; its genomes are marked NOT_DEREPLICATED."
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
log "message" "Dereplication: $(count_lines "${WORK_DIR}/single_leaves.txt") leaf directories with a single genome; ${N_MULTI} with several genomes (${DEREP_METHOD})."

# 2. Leaves with several genomes: skDER or CiDDER, one leaf at a time with every CPU
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
