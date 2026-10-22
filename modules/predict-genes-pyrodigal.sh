#!/usr/bin/env bash
#
# predict-genes-pyrodigal.sh - Prokfinder module 5: gene prediction of the genomes (Pyrodigal)
#
# The genes of each genome that passed the structural evaluation are predicted once, with Pyrodigal (the Python
# port of Prodigal), and the predicted proteins are used by every next step: CheckM2 (--genes), BUSCO (proteins
# mode), GTDB-Tk (identify --genes) and CiDDER (--proteomes). The work is made by predict-genes-pyrodigal.py:
#   - single mode, with the translation table chosen as GTDB-Tk does: table 4 when its coding density is more than
#     5 points higher than with table 11 and above 70%, otherwise table 11
#   - meta mode only for genomes shorter than 100 kbp (too short to train Prodigal). A MAG is the genome of a single
#     organism, so it is predicted in single mode too, as GTDB-Tk does
#
# Re-execution: the result of each genome is kept in a cache table with the modification time of the genome file, so
# only new or changed genomes are predicted again.
#
# Usage: predict-genes-pyrodigal.sh -l <accession_list> -i <clean_genomes_dir> -o <genes_dir> -t <table> [-p <threads>]
#
# Output:
#   <genes_dir>/<ACCESSION>.faa / .ffn / .gff   proteins, gene sequences and coordinates
#   <genes_dir>/gene-prediction-cache.tsv
#   <table>: ASSEMBLY_ACCESSION GENE_CALLING_MODE TRANSLATION_TABLE CODING_DENSITY_TABLE11 CODING_DENSITY_TABLE4
#            PREDICTED_GENES PROTEIN_FILE GENE_PREDICTION_STATUS (PREDICTED | FAILED)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -i <clean_genomes_dir> -o <genes_dir> -t <table> [-p <threads>]

  -l <file>         Accessions whose genes are predicted (one per line)
  -i <dir>          Folder with the clean genomes (<ACCESSION>.fna or .fna.gz)
  -o <dir>          Folder of the predicted genes (.faa, .ffn, .gff)
  -t <file>         Result table (one row per genome)
  -p <n>            Genomes predicted at the same time (default: 4)

Needs Python 3 with the package pyrodigal (on the UFV cluster, the Python of the skder module has it).
EOF
    exit "${1:-1}"
}

ACC_LIST=""
CLEAN_DIR=""
GENES_DIR=""
TABLE=""
THREADS=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -i) CLEAN_DIR="${2%/}"; shift 2 ;;
        -o) GENES_DIR="${2%/}"; shift 2 ;;
        -t) TABLE="$2"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${ACC_LIST}" || -z "${CLEAN_DIR}" || -z "${GENES_DIR}" || -z "${TABLE}" ]]; then
    log "error" "Options -l, -i, -o and -t are required."
    usage
fi
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"
[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || die "The number of threads must be a positive integer (got '${THREADS}')."

mkdir -p "${GENES_DIR}"
setup_module_log "${GENES_DIR}/predict-genes-pyrodigal.log"

PYTHON="$(pf_python_with_pyrodigal)"
[[ -n "${PYTHON}" ]] || die "No Python with the package pyrodigal was found (install it: pip install pyrodigal, or conda install -c bioconda pyrodigal; on the UFV cluster: module load skder)."

CACHE="${GENES_DIR}/gene-prediction-cache.tsv"
# Gene calling parameters recorded in the cache ("meta_for_mags=0", written by older versions, is the same prediction)
PARAMS="default"
LEGACY_PARAMS="meta_for_mags=0"
HEADER='ASSEMBLY_ACCESSION\tFILE_MTIME\tPARAMETERS\tGENE_CALLING_MODE\tTRANSLATION_TABLE\tCODING_DENSITY_TABLE11\tCODING_DENSITY_TABLE4\tPREDICTED_GENES\tPROTEIN_FILE\n'
[[ -s "${CACHE}" ]] || printf "${HEADER}" > "${CACHE}"

make_temp_dir WORK_DIR
GENOMES="${WORK_DIR}/genomes.tsv"   # <accession> <TAB> <file> <TAB> <mtime>

# 1. Genome files of the requested accessions
while IFS= read -r ACC; do
    FILE="$(genome_file_of "${CLEAN_DIR}" "${ACC}")"
    if [[ -n "${FILE}" ]]; then
        printf '%s\t%s\t%s\n' "${ACC}" "$(absolute_path "${FILE}")" "$(file_mtime "${FILE}")"
    fi
done < <(tr -d '\r' < "${ACC_LIST}" | awk 'NF { print $1 }' | sort -u) > "${GENOMES}"

N_GENOMES="$(count_lines "${GENOMES}")"
[[ "${N_GENOMES}" -gt 0 ]] || die "No genome of the requested accessions was found in ${CLEAN_DIR}."

# 2. Genomes to predict: not in the cache with the same mtime and options, or whose .faa file is missing
awk -F'\t' -v params="${PARAMS}" -v legacy="${LEGACY_PARAMS}" -v dir="${GENES_DIR}" '
    NR == FNR { if (FNR > 1 && ($3 == params || $3 == legacy)) done[$1 SUBSEP $2] = 1; next }
    { faa = dir "/" $1 ".faa" }
    !(($1 SUBSEP $3) in done) || (system("test -s \"" faa "\"") != 0)
' "${CACHE}" "${GENOMES}" > "${WORK_DIR}/pending.tsv"
N_PENDING="$(count_lines "${WORK_DIR}/pending.tsv")"
log "message" "Gene prediction (Pyrodigal): ${N_GENOMES} genome(s); $((N_GENOMES - N_PENDING)) reused; ${N_PENDING} to predict (${THREADS} at a time)."
log "message" "  Mode: single, translation table 11 or 4 (GTDB-Tk heuristic); meta for genomes < 100 kbp."

if [[ "${N_PENDING}" -gt 0 ]]; then
    log "message" "  Python: ${PYTHON} ($("${PYTHON}" -c 'import pyrodigal; print("pyrodigal", pyrodigal.__version__)'))"
    cut -f1,2 "${WORK_DIR}/pending.tsv" > "${WORK_DIR}/pending_genomes.tsv"
    PY_ARGS=(--genomes "${WORK_DIR}/pending_genomes.tsv" --out-dir "${GENES_DIR}" --table "${WORK_DIR}/predicted.tsv"
             --threads "${THREADS}")
    "${PYTHON}" "${SCRIPT_DIR}/predict-genes-pyrodigal.py" "${PY_ARGS[@]}" || die "predict-genes-pyrodigal.py failed."

    # New cache rows: <accession> <mtime> <parameters> <results...>
    awk -F'\t' -v OFS='\t' -v params="${PARAMS}" '
        NR == FNR { mtime[$1] = $3; next }
        FNR > 1 { acc = $1; $1 = ""; sub(/^\t/, ""); print acc, mtime[acc], params, $0 }
    ' "${WORK_DIR}/pending.tsv" "${WORK_DIR}/predicted.tsv" > "${WORK_DIR}/new_rows.tsv"
    update_cache "${CACHE}" "${WORK_DIR}/new_rows.tsv"
fi

# 3. Result table (genomes without prediction: FAILED)
awk -F'\t' -v OFS='\t' -v params="${PARAMS}" -v legacy="${LEGACY_PARAMS}" '
    NR == FNR { if (FNR > 1 && ($3 == params || $3 == legacy)) row[$1 SUBSEP $2] = $4 OFS $5 OFS $6 OFS $7 OFS $8 OFS $9; next }
    FNR == 1 { print "ASSEMBLY_ACCESSION", "GENE_CALLING_MODE", "TRANSLATION_TABLE", "CODING_DENSITY_TABLE11",
               "CODING_DENSITY_TABLE4", "PREDICTED_GENES", "PROTEIN_FILE", "GENE_PREDICTION_STATUS"; next }
    {
        key = $1 SUBSEP $3
        if (key in row) print $1, row[key], "PREDICTED"
        else print $1, "NULL", "NULL", "NULL", "NULL", "NULL", "NULL", "FAILED"
    }
' "${CACHE}" <(printf 'header\n'; cat "${GENOMES}") > "${TABLE}"

N_OK="$(awk -F'\t' 'NR > 1 && $8 == "PREDICTED"' "${TABLE}" | wc -l)"
N_TABLE4="$(awk -F'\t' 'NR > 1 && $8 == "PREDICTED" && $3 == "4"' "${TABLE}" | wc -l)"
N_META="$(awk -F'\t' 'NR > 1 && $2 == "meta"' "${TABLE}" | wc -l)"
log "message" "Genes predicted for ${N_OK} of ${N_GENOMES} genome(s) (${N_TABLE4} with translation table 4; ${N_META} in meta mode). Table: ${TABLE}"
if [[ "${N_OK}" -lt "${N_GENOMES}" ]]; then
    log "warning" "$((N_GENOMES - N_OK)) genome(s) without predicted genes (GENE_PREDICTION_STATUS = FAILED) are left out of the next steps."
fi
