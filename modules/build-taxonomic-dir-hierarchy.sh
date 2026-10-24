#!/usr/bin/env bash
#
# build-taxonomic-dir-hierarchy.sh - Prokfinder module 8: taxonomic directory hierarchy
#
# Builds the directory tree Domain/Phylum/Class/Order/Family/Genus/Species from the GTDB taxonomy of the genomes.
# The leaf directory of a genome is its lowest known rank (species or, when the species is unknown, genus,
# family, ...). The genome files are not copied or moved: each leaf directory gets a list with the paths of its
# genomes (genomes-files-symlink.txt), which the dereplication step reads.
# Genomes without GTDB taxonomy (NON-CLASSIFIED-GENOMES-BY-GTDB.txt) are not placed in the hierarchy.
#
# Re-execution: the lists are rebuilt on every run (it is fast), so genomes that left the selection also leave
# the hierarchy. The dereplication results are not kept here but in OUTPUT/GENOME-DEREPLICATION (and the
# REPRESENTATIVE_GENOME column only in genome-metadata.tsv); the GENOME-DEREPLICATION and GENOME-METADATA folders
# written in the leaves by older versions are removed.
#
# Usage: build-taxonomic-dir-hierarchy.sh -m <genome_metadata.tsv> -l <accession_list> -o <hierarchy_dir>
#   <genome_metadata.tsv> must have the columns GTDB_DOMAIN ... GTDB_SPECIES and CLEAN_GENOME_FILE
#
# Output:
#   <hierarchy_dir>/<Domain>/.../<Species>/genomes-files-symlink.txt   (ASSEMBLY_ACCESSION <TAB> genome path)
#   <hierarchy_dir>/genome-dir-hierarchy.tsv   ASSEMBLY_ACCESSION  TAXONOMIC_LEAF_RANK  TAXONOMIC_LEAF_DIR
#   <hierarchy_dir>/leaf-directories.tsv       TAXONOMIC_LEAF_DIR  TAXONOMIC_LEAF_RANK  NUMBER_OF_GENOMES

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -m <genome_metadata.tsv> -l <accession_list> -o <hierarchy_dir>

  -m <file>   Genome metadata table with the GTDB_<RANK> columns and CLEAN_GENOME_FILE
  -l <file>   Accessions placed in the hierarchy (one per line)
  -o <dir>    Root of the directory hierarchy
EOF
    exit "${1:-1}"
}

METADATA=""
ACC_LIST=""
HIERARCHY_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m) METADATA="$2"; shift 2 ;;
        -l) ACC_LIST="$2"; shift 2 ;;
        -o) HIERARCHY_DIR="${2%/}"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${METADATA}" || -z "${ACC_LIST}" || -z "${HIERARCHY_DIR}" ]]; then
    log "error" "Options -m, -l and -o are required."
    usage
fi
[[ -f "${METADATA}" ]] || die "Metadata table not found: ${METADATA}"
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"

mkdir -p "${HIERARCHY_DIR}"
setup_module_log "${HIERARCHY_DIR}/build-taxonomic-dir-hierarchy.log"
HIERARCHY_DIR="$(absolute_path "${HIERARCHY_DIR}")"

GENOME_TABLE="${HIERARCHY_DIR}/genome-dir-hierarchy.tsv"
LEAF_TABLE="${HIERARCHY_DIR}/leaf-directories.tsv"
LIST_NAME="genomes-files-symlink.txt"

RANK_COLS=""
for RANK in DOMAIN PHYLUM CLASS ORDER FAMILY GENUS SPECIES; do
    RANK_COLS="${RANK_COLS:+${RANK_COLS},}$(require_column "${METADATA}" "GTDB_${RANK}")"
done
C_FILE="$(require_column "${METADATA}" CLEAN_GENOME_FILE)"

make_temp_dir WORK_DIR

# 1. Leaf directory of each genome: the ranks are followed from the domain down and the path stops at the first
#    rank without a name. Names are made safe for directory names ("/" and spaces -> "_").
awk -F'\t' -v OFS='\t' -v acc_file="${ACC_LIST}" -v rank_cols="${RANK_COLS}" -v c_file="${C_FILE}" -v root="${HIERARCHY_DIR}" '
    BEGIN {
        while ((getline a < acc_file) > 0) { gsub(/\r/, "", a); if (a != "") want[a] = 1 }
        close(acc_file)
        nr = split(rank_cols, rc, ",")
        split("DOMAIN PHYLUM CLASS ORDER FAMILY GENUS SPECIES", rank_name, " ")
    }
    NR == 1 { next }
    ($1 in want) {
        path = ""; leaf_rank = ""
        for (k = 1; k <= nr; k++) {
            v = $(rc[k])
            if (v == "" || v == "NULL") break
            gsub(/[\/[:space:]]+/, "_", v)
            path = path "/" v
            leaf_rank = rank_name[k]
        }
        if (path == "" || $c_file == "" || $c_file == "NULL") { print $1 > "/dev/stderr"; next }
        print $1, leaf_rank, root path, $c_file
    }
' "${METADATA}" 2> "${WORK_DIR}/not_placed.txt" | LC_ALL=C sort -t $'\t' -k3,3 -k1,1 > "${WORK_DIR}/placement.tsv"

N_PLACED="$(count_lines "${WORK_DIR}/placement.tsv")"
if [[ -s "${WORK_DIR}/not_placed.txt" ]]; then
    log "warning" "$(count_lines "${WORK_DIR}/not_placed.txt") genome(s) have no GTDB taxonomy or no clean genome file and were not placed in the hierarchy."
fi
log "message" "Placing ${N_PLACED} genome(s) in the taxonomic directory hierarchy: ${HIERARCHY_DIR}"

# 2. Remove the lists of the previous run, then create the leaf directories and write their lists
find "${HIERARCHY_DIR}" -mindepth 2 -type f -name "${LIST_NAME}" -delete

cut -f3 "${WORK_DIR}/placement.tsv" | uniq | xargs -r -d '\n' mkdir -p

awk -F'\t' -v list_name="${LIST_NAME}" '
    $3 != current { if (current != "") close(out); current = $3; out = current "/" list_name }
    { print $1 "\t" $4 > out }
' "${WORK_DIR}/placement.tsv"

# 3. Clean-up: folders written in the leaves by older versions (dereplication results and per-leaf metadata), then
#    empty directories
find "${HIERARCHY_DIR}" -mindepth 2 -type d \( -name "GENOME-DEREPLICATION" -o -name "GENOME-METADATA" \) -prune \
    -exec rm -rf {} +
find "${HIERARCHY_DIR}" -mindepth 1 -type d -empty -delete

# 4. Tables
{
    printf 'ASSEMBLY_ACCESSION\tTAXONOMIC_LEAF_RANK\tTAXONOMIC_LEAF_DIR\n'
    cut -f1-3 "${WORK_DIR}/placement.tsv" | LC_ALL=C sort -t $'\t' -k1,1
} > "${GENOME_TABLE}"

{
    printf 'TAXONOMIC_LEAF_DIR\tTAXONOMIC_LEAF_RANK\tNUMBER_OF_GENOMES\n'
    awk -F'\t' -v OFS='\t' '{ n[$3]++; r[$3] = $2 } END { for (d in n) print d, r[d], n[d] }' "${WORK_DIR}/placement.tsv" | LC_ALL=C sort
} > "${LEAF_TABLE}"

log "message" "Hierarchy complete: $(count_records "${LEAF_TABLE}") leaf directories (by rank: $(awk -F'\t' 'NR > 1 { c[$2]++ } END { for (r in c) printf "%s=%d ", r, c[r] }' "${LEAF_TABLE}"))"
log "message" "Tables: ${GENOME_TABLE} and ${LEAF_TABLE}"
