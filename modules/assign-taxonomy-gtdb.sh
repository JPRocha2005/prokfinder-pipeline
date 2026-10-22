#!/usr/bin/env bash
#
# assign-taxonomy-gtdb.sh - Prokfinder module 7: GTDB taxonomy of the genomes
#
# Each genome gets the GTDB taxonomy from the first of these methods that classifies it (GTDB_CLASSIFICATION_METHOD):
#   GTDB_TABLE      the genome is in the GTDB taxonomy tables of the reference databases. GTDB lists each genome under
#                   only one of its accessions, so the lookup tries, in this order: the same accession; the paired
#                   GenBank/RefSeq accession (GCA_x.1 <-> GCF_x.1); another version (GCA_x.*)
#   GTDBTK_ANI      gtdbtk ani_rep (skani, on the genome sequences): the genome meets the GTDB species circumscription
#                   criteria with its closest GTDB representative genome (ANI >= the ANI radius of that species,
#                   usually 95%, and aligned fraction >= 0.5) and gets the taxonomy of that species. It needs GTDB-Tk
#                   and its reference package (GTDBTK_DATA_PATH), but little RAM; without them this step is skipped
#   GTDBTK_TREE     only with --use-gtdbtk-tree: the genomes left are placed in the GTDB reference tree, from the
#                   proteins predicted by Pyrodigal (gtdbtk identify --genes -> align -> classify --genes; with --genes,
#                   GTDB-Tk skips its own ANI step, which GTDBTK_ANI replaces). The placement needs a lot of RAM
#   NOT_CLASSIFIED  none of them; listed in NON-CLASSIFIED-GENOMES-BY-GTDB.txt and left out of the next steps
#
# Re-execution: the ANI and tree results of each genome are kept in cache tables (GTDBTK-RESULTS/), so GTDB-Tk runs
# only once per genome and reference package. The same command run again with --use-gtdbtk-tree only places in the
# tree the genomes that are still NOT_CLASSIFIED.
#
# Usage: assign-taxonomy-gtdb.sh -l <accession_list> -g <gtdb_dir> -o <taxonomy_dir> [-i <clean_genomes_dir>]
#                                [-a <genes_dir>] [-p <threads>] [--memory <GB>] [--use-gtdbtk-tree]
#
# Output inside <taxonomy_dir>:
#   genome-taxonomy.tsv                 ASSEMBLY_ACCESSION GTDB_DOMAIN ... GTDB_SPECIES GTDB_CLASSIFICATION_METHOD
#                                       GTDB_TAXONOMY_SOURCE and the status of the genome in GTDB (added to
#                                       genome-metadata.tsv): GTDB_SPECIES_REPRESENTATIVE (yes/no),
#                                       GTDB_REPRESENTATIVE_GENOME (representative of its species), GTDB_TYPE_MATERIAL,
#                                       GTDB_TYPE_SPECIES_OF_GENUS. A genome classified by ANI is not a representative
#                                       (no) and its representative is the closest GTDB representative genome
#   taxonomy-classification.tsv         how each genome was classified: method, source (GTDB release or GTDB-Tk
#                                       version and reference package), accession found in the GTDB table, closest
#                                       GTDB representative genome, ANI, aligned fraction, the thresholds used, method
#                                       of the GTDB-Tk tree placement and a note (why a genome was not classified)
#   NON-CLASSIFIED-GENOMES-BY-GTDB.txt  accessions without GTDB taxonomy
#   GTDBTK-RESULTS/                     GTDB-Tk runs and caches (gtdbtk-ani-cache.tsv, gtdbtk-tree-cache.tsv)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -l <accession_list> -g <gtdb_dir> -o <taxonomy_dir> [-i <clean_genomes_dir>] [-a <genes_dir>]
       [-p <threads>] [--memory <GB>] [--use-gtdbtk-tree]

  -l <file>           Accessions to classify (one per line)
  -g <dir>            ~/prokfinder-refdb/GTDB
  -o <dir>            Output folder
  -i <dir>            Folder with the clean genomes (ANI to the GTDB representatives with GTDB-Tk)
  -a <dir>            Folder with the predicted proteins, <ACCESSION>.faa (required with --use-gtdbtk-tree)
  -p <n>              CPUs used by GTDB-Tk (default: 4)
  --memory <GB>       RAM available (below 200 GB, the tree placement writes to disk to use less RAM)
  --use-gtdbtk-tree   Place in the GTDB tree the genomes that are neither in the GTDB tables nor classified by ANI

GTDBTK_DATA_PATH must point to the GTDB-Tk reference package (the gtdbtk module of the cluster sets it).
EOF
    exit "${1:-1}"
}

ACC_LIST=""
GTDB_DIR=""
TAXONOMY_DIR=""
CLEAN_DIR=""
GENES_DIR=""
USE_TREE=0
THREADS=4
MEMORY_GB=""
TREE_SCRATCH_BELOW_GB=200   # less RAM than this: classify --scratch_dir (pplacer writes to disk; slower)

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l) ACC_LIST="$2"; shift 2 ;;
        -g) GTDB_DIR="${2%/}"; shift 2 ;;
        -o) TAXONOMY_DIR="${2%/}"; shift 2 ;;
        -i) CLEAN_DIR="${2%/}"; shift 2 ;;
        -a) GENES_DIR="${2%/}"; shift 2 ;;
        -p) THREADS="$2"; shift 2 ;;
        --memory) MEMORY_GB="$2"; shift 2 ;;
        --use-gtdbtk-tree) USE_TREE=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${ACC_LIST}" || -z "${GTDB_DIR}" || -z "${TAXONOMY_DIR}" ]]; then
    log "error" "Options -l, -g and -o are required."
    usage
fi
[[ -f "${ACC_LIST}" ]] || die "Accession list not found: ${ACC_LIST}"
if [[ "${USE_TREE}" -eq 1 && ( -z "${CLEAN_DIR}" || -z "${GENES_DIR}" ) ]]; then
    die "Options -i (clean genomes folder) and -a (predicted proteins folder) are required with --use-gtdbtk-tree."
fi
[[ -z "${MEMORY_GB}" || "${MEMORY_GB}" =~ ^[0-9]+$ ]] || die "--memory must be a number of GB (got '${MEMORY_GB}')."

mkdir -p "${TAXONOMY_DIR}"
setup_module_log "${TAXONOMY_DIR}/assign-taxonomy-gtdb.log"

TAXONOMY_TABLE="${TAXONOMY_DIR}/genome-taxonomy.tsv"
DETAILS_TABLE="${TAXONOMY_DIR}/taxonomy-classification.tsv"
NON_CLASSIFIED_LIST="${TAXONOMY_DIR}/NON-CLASSIFIED-GENOMES-BY-GTDB.txt"
GTDBTK_DIR="${TAXONOMY_DIR}/GTDBTK-RESULTS"
ANI_CACHE="${GTDBTK_DIR}/gtdbtk-ani-cache.tsv"
TREE_CACHE="${GTDBTK_DIR}/gtdbtk-tree-cache.tsv"
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

############ 1. GTDB TABLES ################
#   Only the rows of the requested genomes are kept in memory.
#   Output: <accession> <TAB> 7 ranks <TAB> <accession found in GTDB> <TAB> EXACT | PAIRED_ACCESSION | OTHER_VERSION
#   <TAB> 4 columns of the status in GTDB (the genomes not found have NULL in every column)

STATUS_COLS="GTDB_SPECIES_REPRESENTATIVE GTDB_REPRESENTATIVE_GENOME GTDB_TYPE_MATERIAL GTDB_TYPE_SPECIES_OF_GENUS"
awk -F'\t' -v OFS='\t' -v acc_file="${WORK_DIR}/accessions.txt" -v ranks="${RANKS} ${STATUS_COLS}" '
    function num(a) { sub(/^GC[AF]_/, "", a); return a }
    function base(a) { a = num(a); sub(/[.][0-9]+$/, "", a); return a }
    BEGIN {
        nr = split(ranks, rank, " ")
        for (k = 1; k <= nr; k++) name[k] = (k <= 7 ? "GTDB_" rank[k] : rank[k])
        while ((getline a < acc_file) > 0) { want[++n] = a; by_exact[a] = 1; by_num[num(a)] = 1; by_base[base(a)] = 1 }
        close(acc_file)
    }
    FNR == 1 {
        delete col
        for (i = 1; i <= NF; i++) for (k = 1; k <= nr; k++) if ($i == name[k]) col[k] = i
        next
    }
    ($1 in by_exact) || (num($1) in by_num) || (base($1) in by_base) {
        ranks7 = ""; status = ""
        for (k = 1; k <= nr; k++) {
            v = ((k in col) ? $(col[k]) : "NULL")
            if (k <= 7) ranks7 = ranks7 OFS v; else status = status OFS v
        }
        line = ranks7 OFS $1 SUBSEP status
        exact[$1] = line; paired[num($1)] = line; other[base($1)] = line
    }
    function out(a, line, match_type,   p) {
        split(line, p, SUBSEP)
        print a p[1], match_type p[2]
    }
    END {
        for (i = 1; i <= n; i++) {
            a = want[i]
            if (a in exact)              out(a, exact[a], "EXACT")
            else if (num(a) in paired)   out(a, paired[num(a)], "PAIRED_ACCESSION")
            else if (base(a) in other)   out(a, other[base(a)], "OTHER_VERSION")
            else {
                line = a
                for (k = 1; k <= nr + 2; k++) line = line OFS "NULL"
                print line
            }
        }
    }
' "${GTDB_TABLES[@]}" > "${WORK_DIR}/table_lookup.tsv"

awk -F'\t' '$10 == "NULL" { print $1 }' "${WORK_DIR}/table_lookup.tsv" > "${WORK_DIR}/not_in_gtdb.txt"
N_NOT_IN_GTDB="$(count_lines "${WORK_DIR}/not_in_gtdb.txt")"
log "message" "Found in the GTDB tables: $(( $(count_lines "${WORK_DIR}/table_lookup.tsv") - N_NOT_IN_GTDB )); not found: ${N_NOT_IN_GTDB}"

############ 2. GTDB-TK ################

GTDBTK_OK=0          # GTDB-Tk and its reference package are available
GTDBTK_VERSION="NULL"
REFPKG="NULL"        # reference package of GTDB-Tk (e.g., r232): the caches are valid for one package
ANI_CACHE_HEADER='ASSEMBLY_ACCESSION\tREFERENCE_PACKAGE\tGTDBTK_VERSION\tANI_REFERENCE_GENOME\tANI\tALIGNED_FRACTION\tANI_THRESHOLD\tAF_THRESHOLD\tSATISFIES_GTDB_CRITERIA\tREFERENCE_TAXONOMY\n'
TREE_CACHE_HEADER='ASSEMBLY_ACCESSION\tREFERENCE_PACKAGE\tGTDBTK_VERSION\tCLASSIFICATION\tTREE_CLASSIFICATION_METHOD\tNOTE\n'

if [[ "${N_NOT_IN_GTDB}" -gt 0 ]]; then
    if command -v gtdbtk >/dev/null 2>&1 && [[ -n "${GTDBTK_DATA_PATH:-}" && -d "${GTDBTK_DATA_PATH}" ]]; then
        GTDBTK_OK=1
        GTDBTK_VERSION="$(run_in_own_env gtdbtk --version 2>&1 | grep -oE '[0-9]+[.][0-9]+[.][0-9]+' | head -n 1 || true)"
        GTDBTK_VERSION="${GTDBTK_VERSION:-unknown}"
        REFPKG="$(basename "${GTDBTK_DATA_PATH%/}")"
        mkdir -p "${GTDBTK_DIR}"
        [[ -s "${ANI_CACHE}" ]] || printf "${ANI_CACHE_HEADER}" > "${ANI_CACHE}"
        [[ -s "${TREE_CACHE}" ]] || printf "${TREE_CACHE_HEADER}" > "${TREE_CACHE}"
        log "message" "GTDB-Tk ${GTDBTK_VERSION}, reference package ${REFPKG} (${GTDBTK_DATA_PATH})"
    elif [[ "${USE_TREE}" -eq 1 ]]; then
        die "--use-gtdbtk-tree needs GTDB-Tk and its reference package (GTDBTK_DATA_PATH); one of them was not found."
    else
        log "warning" "GTDB-Tk or its reference package (GTDBTK_DATA_PATH) was not found: the ${N_NOT_IN_GTDB} genome(s) not in the GTDB tables are not classified by ANI."
    fi
fi

# Runs a GTDB-Tk command (its output goes to the console log of the run); stops on failure
# Usage: run_gtdbtk <console_log> gtdbtk <step> [arguments]
run_gtdbtk() {
    local CONSOLE_LOG="$1"
    shift
    echo "### $*" >> "${CONSOLE_LOG}"
    if ! run_in_own_env "$@" >> "${CONSOLE_LOG}" 2>&1; then
        log "error" "GTDB-Tk failed ($2). Last lines of ${CONSOLE_LOG}:"
        tail -n 20 "${CONSOLE_LOG}" >&2
        exit 1
    fi
}

# Genomes of <list> without a row of the current reference package in <cache>
# Usage: not_in_cache <cache> <list>
not_in_cache() {
    awk -F'\t' -v refpkg="${REFPKG}" 'NR == FNR { if (FNR > 1 && $2 == refpkg) done[$1] = 1; next } !($1 in done)' "$1" "$2"
}

# 2a. ANI to the GTDB representative genomes (gtdbtk ani_rep), for the genomes not in the GTDB tables
if [[ "${GTDBTK_OK}" -eq 1 ]]; then
    not_in_cache "${ANI_CACHE}" "${WORK_DIR}/not_in_gtdb.txt" > "${WORK_DIR}/ani_pending.txt"
    # Batch file of GTDB-Tk: <genome file> <TAB> <genome id>
    : > "${WORK_DIR}/ani_batch.tsv"
    while IFS= read -r ACC; do
        FILE="$( [[ -z "${CLEAN_DIR}" ]] || genome_file_of "${CLEAN_DIR}" "${ACC}")"
        [[ -z "${FILE}" ]] || printf '%s\t%s\n' "$(absolute_path "${FILE}")" "${ACC}" >> "${WORK_DIR}/ani_batch.tsv"
    done < "${WORK_DIR}/ani_pending.txt"
    N_ANI_RUN="$(count_lines "${WORK_DIR}/ani_batch.tsv")"
    log "message" "GTDB-Tk ani_rep: $(( N_NOT_IN_GTDB - $(count_lines "${WORK_DIR}/ani_pending.txt") )) genome(s) reused from cache; ${N_ANI_RUN} to compare with the GTDB representative genomes."

    if [[ "${N_ANI_RUN}" -gt 0 ]]; then
        RUN_DIR="${GTDBTK_DIR}/ani_rep_$(date +%Y%m%d_%H%M%S)"
        mkdir -p "${RUN_DIR}/tmp"
        log "message" "  gtdbtk ani_rep: ANI of ${N_ANI_RUN} genome(s) to the GTDB representative genomes (${THREADS} CPUs)..."
        run_gtdbtk "${RUN_DIR}/gtdbtk-console.log" gtdbtk ani_rep --batchfile "${WORK_DIR}/ani_batch.tsv" \
            --out_dir "${RUN_DIR}" --cpus "${THREADS}" --tmpdir "${RUN_DIR}/tmp"
        rm -rf "${RUN_DIR}/tmp"
        ANI_CLOSEST="$(find "${RUN_DIR}" -type f -name '*ani_closest.tsv' | head -n 1)"
        [[ -n "${ANI_CLOSEST}" ]] || die "GTDB-Tk ani_rep wrote no *ani_closest.tsv file (see ${RUN_DIR}/gtdbtk-console.log)."

        # ani_closest.tsv: user_genome, reference_genome, skani_ani, skani_af, reference_taxonomy,
        # satisfies_gtdb_circumscription_criteria ("no result" when no representative has AF >= 0.5). The ANI
        # threshold is the ANI radius of the species of the representative (radii/gtdb_radii.tsv of the reference
        # package: species, representative genome, radius); the AF threshold is 0.5 (GTDB-Tk AF_THRESHOLD)
        RADII="${GTDBTK_DATA_PATH%/}/radii/gtdb_radii.tsv"
        [[ -f "${RADII}" ]] || RADII="/dev/null"
        awk -F'\t' -v OFS='\t' -v refpkg="${REFPKG}" -v version="${GTDBTK_VERSION}" -v batch="${WORK_DIR}/ani_batch.tsv" \
            -v radii="${RADII}" '
            function canonical(g) { sub(/^(RS_|GB_)/, "", g); return g }
            function value(x) { return (x == "" || x == "no result" || x == "N/A") ? "NULL" : x }
            BEGIN {
                while ((getline line < radii) > 0) { split(line, f, "\t"); if (f[3] ~ /^[0-9.]+$/) radius[canonical(f[2])] = f[3] + 0 }
                close(radii)
            }
            FNR == 1 {
                for (i = 1; i <= NF; i++) {
                    h = tolower($i)
                    if (h == "user_genome") cu = i
                    else if (h == "reference_genome") cr = i
                    else if (h ~ /satisf/) cs = i
                    else if (h ~ /taxonomy/) ct = i
                    else if (h ~ /(^|_)ani$/) ca = i
                    else if (h ~ /(^|_)af$/) cf = i
                }
                if (!cu || !ct || !cs) { print "ani_closest.tsv: columns user_genome, *taxonomy* and satisfies_* are required" > "/dev/stderr"; exit 2 }
                next
            }
            {
                if (!match($cu, /GC[AF]_[0-9]+[.][0-9]+/)) next
                acc = substr($cu, RSTART, RLENGTH)
                if (acc in seen) next
                seen[acc] = 1
                ref = (cr ? value($cr) : "NULL")
                thr = (ref != "NULL" && (canonical(ref) in radius)) ? sprintf("%.2f", radius[canonical(ref)]) : "NULL"
                ok = tolower($cs) == "true" ? "TRUE" : (value($cs) == "NULL" ? "NO_HIT" : "FALSE")
                print acc, refpkg, version, ref, (ca ? value($ca) : "NULL"), (cf ? value($cf) : "NULL"), thr, "0.5", ok, value($ct)
            }
            END {
                # genomes of the batch without a row: no representative genome close enough
                while ((getline line < batch) > 0) {
                    split(line, f, "\t")
                    if (!(f[2] in seen)) print f[2], refpkg, version, "NULL", "NULL", "NULL", "NULL", "0.5", "NO_HIT", "NULL"
                }
            }
        ' "${ANI_CLOSEST}" > "${WORK_DIR}/ani_new_rows.tsv" || die "Could not read ${ANI_CLOSEST} (its columns may have changed)."
        update_cache "${ANI_CACHE}" "${WORK_DIR}/ani_new_rows.tsv"
        log "message" "  Classified by ANI to a GTDB representative genome: $(awk -F'\t' '$9 == "TRUE"' "${WORK_DIR}/ani_new_rows.tsv" | wc -l) of ${N_ANI_RUN}."
    fi
fi

# 2b. Placement in the GTDB tree (--use-gtdbtk-tree) of the genomes not classified by ANI
if [[ "${GTDBTK_OK}" -eq 1 ]]; then
    awk -F'\t' -v refpkg="${REFPKG}" 'NR == FNR { if (FNR > 1 && $2 == refpkg && $9 == "TRUE") ani[$1] = 1; next } !($1 in ani)' \
        "${ANI_CACHE}" "${WORK_DIR}/not_in_gtdb.txt" > "${WORK_DIR}/tree_candidates.txt"
    N_TREE_CANDIDATES="$(count_lines "${WORK_DIR}/tree_candidates.txt")"
    not_in_cache "${TREE_CACHE}" "${WORK_DIR}/tree_candidates.txt" > "${WORK_DIR}/tree_pending.txt"
    N_TREE_PENDING="$(count_lines "${WORK_DIR}/tree_pending.txt")"

    if [[ "${USE_TREE}" -eq 0 ]]; then
        if [[ "${N_TREE_PENDING}" -gt 0 ]]; then
            log "message" "${N_TREE_PENDING} genome(s) are neither in the GTDB tables nor classified by ANI: they can be placed in the GTDB tree by running the same command again with download-wf --use-gtdbtk-tree (needs a lot of RAM)."
        fi
    elif [[ "${N_TREE_CANDIDATES}" -gt 0 ]]; then
        # Batch file of GTDB-Tk with the predicted proteins: <.faa> <TAB> <genome id>
        : > "${WORK_DIR}/tree_batch.tsv"
        while IFS= read -r ACC; do
            [[ ! -s "${GENES_DIR}/${ACC}.faa" ]] || printf '%s\t%s\n' "$(absolute_path "${GENES_DIR}/${ACC}.faa")" "${ACC}" >> "${WORK_DIR}/tree_batch.tsv"
        done < "${WORK_DIR}/tree_pending.txt"
        N_TREE_RUN="$(count_lines "${WORK_DIR}/tree_batch.tsv")"
        log "message" "GTDB-Tk tree placement (--use-gtdbtk-tree): $(( N_TREE_CANDIDATES - N_TREE_PENDING )) genome(s) reused from cache; ${N_TREE_RUN} to place."

        if [[ "${N_TREE_RUN}" -gt 0 ]]; then
            RUN_DIR="${GTDBTK_DIR}/tree_$(date +%Y%m%d_%H%M%S)"
            CONSOLE_LOG="${RUN_DIR}/gtdbtk-console.log"
            mkdir -p "${RUN_DIR}"
            CLASSIFY_ARGS=()
            if [[ -n "${MEMORY_GB}" ]] && (( MEMORY_GB < TREE_SCRATCH_BELOW_GB )); then
                mkdir -p "${RUN_DIR}/scratch"
                CLASSIFY_ARGS=(--scratch_dir "${RUN_DIR}/scratch")
                log "message" "  ${MEMORY_GB} GB of RAM (less than ${TREE_SCRATCH_BELOW_GB} GB): pplacer writes to disk to use less RAM (--scratch_dir; slower)."
            fi
            log "message" "  gtdbtk identify --genes / align / classify --genes for ${N_TREE_RUN} genome(s) (${THREADS} CPUs)..."
            run_gtdbtk "${CONSOLE_LOG}" gtdbtk identify --batchfile "${WORK_DIR}/tree_batch.tsv" --genes \
                --out_dir "${RUN_DIR}/identify" --cpus "${THREADS}"
            run_gtdbtk "${CONSOLE_LOG}" gtdbtk align --identify_dir "${RUN_DIR}/identify" --out_dir "${RUN_DIR}/align" \
                --cpus "${THREADS}"
            run_gtdbtk "${CONSOLE_LOG}" gtdbtk classify --batchfile "${WORK_DIR}/tree_batch.tsv" --genes \
                --align_dir "${RUN_DIR}/align" --out_dir "${RUN_DIR}/classify" \
                --cpus "${THREADS}" --pplacer_cpus "$(( THREADS > 1 ? THREADS / 2 : 1 ))" \
                ${CLASSIFY_ARGS[@]+"${CLASSIFY_ARGS[@]}"}
            rm -rf "${RUN_DIR}/scratch"

            # gtdbtk.*.summary.tsv: user_genome, classification, classification_method, note...; the genomes without a
            # row (e.g., too few marker genes) are cached as NOT_PLACED, so they are not run again
            find "${RUN_DIR}/classify" -maxdepth 2 -type f -name 'gtdbtk.*.summary.tsv' -print0 \
                | xargs -0 -r awk -F'\t' -v OFS='\t' -v refpkg="${REFPKG}" -v version="${GTDBTK_VERSION}" '
                    FNR == 1 { cg = cc = cm = cn = 0; for (i = 1; i <= NF; i++) { if ($i == "user_genome") cg = i; if ($i == "classification") cc = i; if ($i == "classification_method") cm = i; if ($i == "note") cn = i }; next }
                    {
                        if (!match($cg, /GC[AF]_[0-9]+[.][0-9]+/)) next
                        acc = substr($cg, RSTART, RLENGTH)
                        if (acc in done) next
                        done[acc] = 1
                        cls = ($cc == "" || $cc ~ /^Unclassified/) ? "NULL" : $cc
                        method = (cm && $cm != "" && $cm != "N/A") ? $cm : "NULL"
                        note = (cn && $cn != "" && $cn != "N/A") ? $cn : "NULL"
                        if (cls == "NULL" && note == "NULL") note = "unclassified by GTDB-Tk"
                        print acc, refpkg, version, cls, method, note
                    }' > "${WORK_DIR}/tree_rows.tsv" || true
            awk -F'\t' -v OFS='\t' -v refpkg="${REFPKG}" -v version="${GTDBTK_VERSION}" \
                'NR == FNR { done[$1] = 1; print; next }
                 !($2 in done) { print $2, refpkg, version, "NULL", "NOT_PLACED", "not placed in the tree by GTDB-Tk (see " console ")" }' \
                console="${CONSOLE_LOG}" "${WORK_DIR}/tree_rows.tsv" "${WORK_DIR}/tree_batch.tsv" > "${WORK_DIR}/tree_new_rows.tsv"
            update_cache "${TREE_CACHE}" "${WORK_DIR}/tree_new_rows.tsv"
            log "message" "  Classified by placement in the GTDB tree: $(awk -F'\t' '$4 != "NULL"' "${WORK_DIR}/tree_new_rows.tsv" | wc -l) of ${N_TREE_RUN}."
        fi
    fi
fi

############ 3. FINAL TABLES ################

[[ -f "${ANI_CACHE}" ]] || printf "${ANI_CACHE_HEADER}" > "${WORK_DIR}/empty_ani.tsv"
[[ -f "${TREE_CACHE}" ]] || printf "${TREE_CACHE_HEADER}" > "${WORK_DIR}/empty_tree.tsv"
ANI_TABLE="${ANI_CACHE}"; [[ -f "${ANI_CACHE}" ]] || ANI_TABLE="${WORK_DIR}/empty_ani.tsv"
TREE_TABLE="${TREE_CACHE}"; [[ -f "${TREE_CACHE}" ]] || TREE_TABLE="${WORK_DIR}/empty_tree.tsv"

awk -F'\t' -v OFS='\t' -v ani_file="${ANI_TABLE}" -v tree_file="${TREE_TABLE}" -v details="${DETAILS_TABLE}" \
    -v release="${GTDB_RELEASE}" -v refpkg="${REFPKG}" -v gtdbtk_ok="${GTDBTK_OK}" -v use_tree="${USE_TREE}" \
    -v clean_dir="${CLEAN_DIR}" '
    # "d__Bacteria;p__...;s__Name" -> 7 ranks without the prefixes, TAB-separated (NULL when empty)
    function ranks_of(lineage,   n, lin, k, v, out) {
        n = split(lineage, lin, ";"); out = ""
        for (k = 1; k <= 7; k++) { v = (k <= n ? lin[k] : ""); sub(/^[a-z]__/, "", v); out = out (k > 1 ? OFS : "") (v == "" ? "NULL" : v) }
        return out
    }
    BEGIN {
        while ((getline line < ani_file) > 0) { split(line, f, "\t"); if (f[2] == refpkg) ani[f[1]] = line }
        close(ani_file)
        while ((getline line < tree_file) > 0) { split(line, f, "\t"); if (f[2] == refpkg) tree[f[1]] = line }
        close(tree_file)
        null7 = "NULL" OFS "NULL" OFS "NULL" OFS "NULL" OFS "NULL" OFS "NULL" OFS "NULL"
        print "ASSEMBLY_ACCESSION", "GTDB_DOMAIN", "GTDB_PHYLUM", "GTDB_CLASS", "GTDB_ORDER", "GTDB_FAMILY", "GTDB_GENUS",
              "GTDB_SPECIES", "GTDB_CLASSIFICATION_METHOD", "GTDB_TAXONOMY_SOURCE", "GTDB_SPECIES_REPRESENTATIVE",
              "GTDB_REPRESENTATIVE_GENOME", "GTDB_TYPE_MATERIAL", "GTDB_TYPE_SPECIES_OF_GENUS"
        print "ASSEMBLY_ACCESSION", "CLASSIFICATION_METHOD", "SOURCE", "GTDB_TABLE_ACCESSION", "GTDB_TABLE_MATCH",
              "ANI_REFERENCE_GENOME", "ANI", "ALIGNED_FRACTION", "ANI_THRESHOLD", "AF_THRESHOLD",
              "TREE_CLASSIFICATION_METHOD", "NOTE" > details
    }
    {
        acc = $1
        a_ref = a_ani = a_af = a_thr = a_afthr = "NULL"; t_method = "NULL"; note = "NULL"; table_acc = table_match = "NULL"
        status = "NULL" OFS "NULL" OFS "NULL" OFS "NULL"
        if (acc in ani) { split(ani[acc], A, "\t"); a_ref = A[4]; a_ani = A[5]; a_af = A[6]; a_thr = A[7]; a_afthr = A[8] }
        if ($9 != "NULL") {
            method = "GTDB_TABLE"; source = "GTDB " release; table_acc = $9; table_match = $10
            ranks = $2 OFS $3 OFS $4 OFS $5 OFS $6 OFS $7 OFS $8
            status = $11 OFS $12 OFS $13 OFS $14
        } else if ((acc in ani) && A[9] == "TRUE") {
            method = "GTDBTK_ANI"; source = "GTDB-Tk " A[3] " ani_rep (" refpkg ")"; ranks = ranks_of(A[10])
            rep = A[4]; sub(/^(RS|GB)_/, "", rep)
            status = "no" OFS rep OFS "NULL" OFS "NULL"
        } else if ((acc in tree) && (split(tree[acc], T, "\t") > 0) && T[4] != "NULL") {
            method = "GTDBTK_TREE"; source = "GTDB-Tk " T[3] " classify (" refpkg ")"; ranks = ranks_of(T[4])
            t_method = T[5]; note = T[6]
        } else {
            method = "NOT_CLASSIFIED"; source = "NULL"; ranks = null7
            if (gtdbtk_ok != 1) note = "not in the GTDB tables; GTDB-Tk or its reference package (GTDBTK_DATA_PATH) not available, so no ANI classification"
            else if (!(acc in ani)) note = "not in the GTDB tables; genome file not found for the GTDB-Tk ANI" (clean_dir == "" ? " (no genome folder given)" : "")
            else {
                if (A[9] == "NO_HIT") note = "not in the GTDB tables; no GTDB representative genome with aligned fraction >= " A[8]
                else note = "not in the GTDB tables; closest GTDB representative " A[4] " below the species thresholds (ANI " A[5] " < " A[7] " or AF " A[6] " < " A[8] ")"
                if (acc in tree) { split(tree[acc], T, "\t"); t_method = T[5]; note = note "; GTDB-Tk tree: " T[6] }
                else if (use_tree != 1) note = note "; use download-wf --use-gtdbtk-tree to place it in the GTDB tree"
                else note = note "; protein file not found for the GTDB-Tk tree"
            }
        }
        print acc, ranks, method, source, status
        print acc, method, source, table_acc, table_match, a_ref, a_ani, a_af, a_thr, a_afthr, t_method, note > details
    }
' "${WORK_DIR}/table_lookup.tsv" > "${TAXONOMY_TABLE}"

awk -F'\t' 'NR > 1 && $9 == "NOT_CLASSIFIED" { print $1 }' "${TAXONOMY_TABLE}" > "${NON_CLASSIFIED_LIST}"

log "message" "GTDB taxonomy per classification method:"
while IFS=$'\t' read -r METHOD COUNT; do
    log "message" "    ${METHOD}: ${COUNT}"
done < <(awk -F'\t' 'NR > 1 { c[$9]++ } END { for (m in c) printf "%s\t%d\n", m, c[m] }' "${TAXONOMY_TABLE}" | LC_ALL=C sort)
log "message" "Taxonomy saved to: ${TAXONOMY_TABLE} (how each genome was classified: ${DETAILS_TABLE})"
log "message" "Genomes without GTDB taxonomy (left out of the directory hierarchy): $(count_lines "${NON_CLASSIFIED_LIST}") -> ${NON_CLASSIFIED_LIST}"
