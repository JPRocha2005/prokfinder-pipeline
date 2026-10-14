#!/usr/bin/env bash
#
# summarize-final-results.sh - Prokfinder module 8: final results and statistics
#
# Reads the genome metadata table (one row per selected genome, with the columns added by every step) and writes:
#   high-quality-selected-genomes.txt           accession + path of the clean genome of every genome that passed the
#                                               quality filters and is representative of its taxon
#   high-quality-selected-genomes-metadata.txt  metadata of those genomes only
#   STATISTICS/TABLES/
#     selection-funnel.tsv                      number of genomes left after each step
#     qc-status.tsv                             genomes per quality result
#     genomes-per-gtdb-<rank>.tsv               selected / high-quality / representative genomes per taxon
#     all-prokaryote-genomes-status.tsv.gz      every genome of the NCBI table with flags for selected,
#                                               high-quality and representative (for comparison studies)
#   STATISTICS/GRAPHS/                          plots made by plot-statistics.R (only if Rscript is installed)
#
# Usage: summarize-final-results.sh -m <genome_metadata.tsv> -r <ncbi_dir> -s <genbank|refseq> -o <final_results_dir>

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -m <genome_metadata.tsv> -r <ncbi_dir> -s <genbank|refseq> -o <final_results_dir>

  -m <file>   Genome metadata table of the selected genomes
  -r <dir>    prokfinder-results/INPUT/REFERENCE-METADATA/NCBI (for the table of all prokaryote genomes)
  -s <src>    NCBI table used in the search: genbank or refseq
  -o <dir>    Output folder (FINAL-RESULTS)
EOF
    exit "${1:-1}"
}

METADATA=""
NCBI_DIR=""
SOURCE="genbank"
FINAL_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m) METADATA="$2"; shift 2 ;;
        -r) NCBI_DIR="${2%/}"; shift 2 ;;
        -s) SOURCE="$2"; shift 2 ;;
        -o) FINAL_DIR="${2%/}"; shift 2 ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${METADATA}" || -z "${NCBI_DIR}" || -z "${FINAL_DIR}" ]]; then
    log "error" "Options -m, -r and -o are required."
    usage
fi
[[ -f "${METADATA}" ]] || die "Metadata table not found: ${METADATA}"

TABLES_DIR="${FINAL_DIR}/STATISTICS/TABLES"
GRAPHS_DIR="${FINAL_DIR}/STATISTICS/GRAPHS"
mkdir -p "${TABLES_DIR}" "${GRAPHS_DIR}"
setup_module_log "${FINAL_DIR}/summarize-final-results.log"

HQ_LIST="${FINAL_DIR}/high-quality-selected-genomes.txt"
HQ_METADATA="${FINAL_DIR}/high-quality-selected-genomes-metadata.txt"

# Prints the index of a column of the metadata table, or 0 if the step that adds it was not run
col() {
    local INDEX
    INDEX="$(get_column_index "${METADATA}" "$1")"
    echo "${INDEX:-0}"
}

C_DOWNLOAD="$(col DOWNLOAD_STATUS)"
C_CONTIG="$(col CONTIG_FILTER_STATUS)"
C_QC="$(col QC_STATUS)"
C_TAX_SOURCE="$(col GTDB_TAXONOMY_SOURCE)"
C_REP="$(col REPRESENTATIVE_GENOME)"
C_CLEAN_FILE="$(col CLEAN_GENOME_FILE)"

# 1. High-quality, representative genomes
awk -F'\t' -v OFS='\t' -v c_qc="${C_QC}" -v c_rep="${C_REP}" -v c_file="${C_CLEAN_FILE}" -v list="${HQ_LIST}" '
    NR == 1 { print; print "ASSEMBLY_ACCESSION", "CLEAN_GENOME_FILE" > list; next }
    c_qc > 0 && c_rep > 0 && $c_qc == "PASSED" && $c_rep == "yes" {
        print
        print $1, (c_file > 0 ? $c_file : "NULL") > list
    }
' "${METADATA}" > "${HQ_METADATA}"
log "message" "High-quality representative genomes: $(count_records "${HQ_METADATA}") -> ${HQ_LIST}"

# 2. Selection funnel
awk -F'\t' -v OFS='\t' -v c_dl="${C_DOWNLOAD}" -v c_ct="${C_CONTIG}" -v c_qc="${C_QC}" -v c_tx="${C_TAX_SOURCE}" -v c_rep="${C_REP}" '
    NR == 1 { next }
    {
        n_sel++
        if (c_dl && $c_dl == "DOWNLOADED") n_dl++
        if (c_ct && $c_ct == "OK") n_ct++
        if (c_qc && $c_qc == "PASSED") n_qc++
        if (c_qc && c_tx && $c_qc == "PASSED" && $c_tx != "NULL") n_tx++
        if (c_rep && $c_rep == "yes") n_rep++
    }
    END {
        print "STEP", "GENOMES"
        print "1_selected_in_reference_db", n_sel + 0
        print "2_downloaded", (c_dl ? n_dl + 0 : "NA")
        print "3_clean_genome_contigs_filtered", (c_ct ? n_ct + 0 : "NA")
        print "4_passed_quality_filters", (c_qc ? n_qc + 0 : "NA")
        print "5_gtdb_taxonomy_assigned", (c_tx ? n_tx + 0 : "NA")
        print "6_representative_after_dereplication", (c_rep ? n_rep + 0 : "NA")
    }
' "${METADATA}" > "${TABLES_DIR}/selection-funnel.tsv"

log "message" "Selection funnel:"
while IFS=$'\t' read -r STEP COUNT; do
    log "message" "    ${STEP}: ${COUNT}"
done < <(tail -n +2 "${TABLES_DIR}/selection-funnel.tsv")

# 3. Quality results
if [[ "${C_QC}" -gt 0 ]]; then
    { printf 'QC_STATUS\tGENOMES\n'
      awk -F'\t' -v c="${C_QC}" 'NR > 1 { n[$c]++ } END { for (s in n) printf "%s\t%d\n", s, n[s] }' "${METADATA}" | LC_ALL=C sort -t $'\t' -k2,2nr
    } > "${TABLES_DIR}/qc-status.tsv"
fi

# 4. Genomes per GTDB taxon (phylum, family, genus and species)
for RANK in PHYLUM FAMILY GENUS SPECIES; do
    C_RANK="$(col "GTDB_${RANK}")"
    if [[ "${C_RANK}" -eq 0 ]]; then
        continue
    fi
    OUT="${TABLES_DIR}/genomes-per-gtdb-${RANK,,}.tsv"
    { printf 'GTDB_%s\tSELECTED_GENOMES_WITH_TAXONOMY\tHIGH_QUALITY_GENOMES\tREPRESENTATIVE_GENOMES\n' "${RANK}"
      awk -F'\t' -v OFS='\t' -v c="${C_RANK}" -v c_qc="${C_QC}" -v c_rep="${C_REP}" '
        NR > 1 && $c != "NULL" {
            n[$c]++
            if (c_qc && $c_qc == "PASSED") hq[$c]++
            if (c_rep && $c_rep == "yes") rep[$c]++
        }
        END { for (t in n) print t, n[t], hq[t] + 0, rep[t] + 0 }' "${METADATA}" | LC_ALL=C sort -t $'\t' -k4,4nr -k2,2nr
    } > "${OUT}"
done

# 5. Every prokaryote genome of the NCBI table, with the selection flags
STATUS_TABLE="${TABLES_DIR}/all-prokaryote-genomes-status.tsv.gz"
NCBI_TABLES=()
for DOMAIN in "${PF_DOMAINS[@]}"; do
    TABLE="$(pf_ncbi_metadata_table "${NCBI_DIR}" "${DOMAIN}" "${SOURCE}")"
    if [[ -f "${TABLE}" ]]; then
        NCBI_TABLES+=("${TABLE}")
    fi
done

if [[ "${#NCBI_TABLES[@]}" -gt 0 ]]; then
    log "message" "Writing the table of all prokaryote genomes (${SOURCE}) with the selection flags..."
    awk -F'\t' -v OFS='\t' -v meta="${METADATA}" -v c_qc="${C_QC}" -v c_rep="${C_REP}" '
        BEGIN {
            while ((getline line < meta) > 0) {
                if (++k == 1) continue
                split(line, f, "\t")
                selected[f[1]] = 1
                qc[f[1]] = (c_qc ? f[c_qc] : "NULL")
                rep[f[1]] = (c_rep ? f[c_rep] : "NULL")
            }
            close(meta)
        }
        FNR == 1 { if (NR == 1) print $0, "SELECTED", "QC_STATUS", "HIGH_QUALITY", "REPRESENTATIVE_GENOME"; next }
        {
            if ($1 in selected) print $0, "yes", qc[$1], (qc[$1] == "PASSED" ? "yes" : "no"), rep[$1]
            else                print $0, "no", "NULL", "NULL", "NULL"
        }
    ' "${NCBI_TABLES[@]}" | gzip -c > "${STATUS_TABLE}"
fi

# 6. Graphs (optional: base R only, no extra package)
if command -v Rscript >/dev/null 2>&1; then
    log "message" "Drawing graphs with R ($(Rscript --version 2>&1 | head -n 1))..."
    if ! Rscript "${SCRIPT_DIR}/plot-statistics.R" "${METADATA}" "${TABLES_DIR}" "${GRAPHS_DIR}"; then
        log "warning" "plot-statistics.R failed; the tables were written, but not every graph."
    fi
else
    log "warning" "Rscript not found: graphs were not drawn (the tables are in ${TABLES_DIR})."
fi

log "message" "Final results saved to: ${FINAL_DIR}"
