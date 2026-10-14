#!/usr/bin/env bash
#
# search-genomes-in-reference-db.sh - Prokfinder module 1: genome search in the reference database
#
# Reads a query file (genome-selection-query.txt) and selects, in the NCBI tables, the genome assemblies that
# satisfy every filter. Filters can use the metadata of NCBI, the phenotypic traits of BacDive and the
# taxonomy/metadata of GTDB:
#
#   Stage 1 (BacDive) -> NCBI TaxIDs of the strains whose traits satisfy the BACDIVE_* filters
#   Stage 2 (GTDB)    -> accessions whose GTDB taxonomy/metadata satisfy the GTDB_* filters
#   Stage 3 (NCBI)    -> assemblies of the NCBI table that satisfy the NCBI filters AND are linked to stage 1/2
#   Annotation        -> BacDive traits of the selected assemblies are added as BACDIVE_* columns
#
# QUERY FILE: one filter per line, "<KEY><TAB><VALUE>". Empty lines and lines starting with "#" are ignored.
# Lines with the same key in the first group are combined with OR; everything else is combined with AND.
#
#   OR within the key:
#     ACCESSION         GCA_/GCF_ accession (version optional). GCA_x.1 also matches its paired GCF_x.1
#     TAXID             NCBI TaxID of the organism or of any rank of its lineage (descendant taxa included)
#     NCBI_TAXON        name of a taxon of the NCBI lineage (any rank, exact match, case-insensitive)
#     GTDB_TAXON        name of a taxon of the GTDB lineage (any rank; "g__Thermus" or "Thermus")
#     ISOLATION_SOURCE  text contained in the isolation source (case-insensitive)
#     COUNTRY           country/place at the start of a word of the geographic location (case-insensitive)
#   Single value:
#     ASSEMBLY_SOURCE   genbank (default) or refseq: NCBI table that is searched
#     GENOME_TYPE       isolate_genome, mag or all (default)
#   AND (each line is one more condition):
#     NCBI_RANGE        <column>=<min>:<max>    numeric column of the NCBI table (either side may be empty)
#     NCBI_MATCH        <column>=<text>         text contained in a column of the NCBI table
#     BACDIVE_RANGE     <column>=<min>:<max>    numeric trait of BacDive (any recorded value inside the range)
#     BACDIVE_CLASS     <trait>=<class>[,<class>...]   trait divided in classes (any of the classes)
#     GTDB_RANGE        <column>=<min>:<max>    numeric column of the GTDB metadata
#     GTDB_MATCH        <column>=<text>         text contained in a column of the GTDB metadata
#
# Usage:
#   search-genomes-in-reference-db.sh -q <query_file> --validate-only
#   search-genomes-in-reference-db.sh -q <query_file> -r <REFERENCE-METADATA dir> -o <selected_metadata.tsv>
#
# Output: <selected_metadata.tsv> (one row per assembly: NCBI columns + BACDIVE_* columns)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -q <query_file> [-r <reference_metadata_dir> -o <output_tsv>] [--validate-only]

  -q <file>          Query file (see the header of this script for the format)
  -r <dir>           prokfinder-results/INPUT/REFERENCE-METADATA
  -o <file>          Output table with the selected genome assemblies
  --validate-only    Only check the query file (no reference table is read)
  -h                 Show this help message
EOF
    exit "${1:-1}"
}

QUERY_FILE=""
REF_DIR=""
OUTPUT_TSV=""
VALIDATE_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -q) QUERY_FILE="$2"; shift 2 ;;
        -r) REF_DIR="${2%/}"; shift 2 ;;
        -o) OUTPUT_TSV="$2"; shift 2 ;;
        --validate-only) VALIDATE_ONLY=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Invalid option: $1"; usage ;;
    esac
done

if [[ -z "${QUERY_FILE}" || ! -f "${QUERY_FILE}" ]]; then
    die "Query file not found: '${QUERY_FILE}'"
fi
if [[ "${VALIDATE_ONLY}" -eq 0 && ( -z "${REF_DIR}" || -z "${OUTPUT_TSV}" ) ]]; then
    log "error" "Options -r and -o are required (unless --validate-only is used)."
    usage
fi

if [[ "${VALIDATE_ONLY}" -eq 0 ]]; then
    setup_module_log "$(dirname "${OUTPUT_TSV}")/search-genomes-in-reference-db.log"
fi

# ---------- Traits of the BacDive table that are divided in classes ----------
# Names are the BacDive table column names; classes are the values written by download-reference-db-bacdive.sh.
# This is the only place to edit to add a trait or a class: the validation and the filtering read from it.
TRAIT_ORDER=(
    cell_shape oxygen_tolerance gram_stain motility flagellum_arrangement colony_color
    pigmentation_production forms_multicellular_complex spore_formation
    energy_source electron_donor carbon_source biosafety_level
)

declare -A TRAIT_CLASSES=(
    [cell_shape]="ovoid spiral coccus curved rod pleomorphic filament other"
    [oxygen_tolerance]="aerobe anaerobe facultative microaerophile"
    [gram_stain]="positive negative"
    [motility]="yes no"
    [flagellum_arrangement]="polar peritrichous gliding"
    [colony_color]="colorless cream yellow white brown pink orange red grey green blue"
    [pigmentation_production]="yes no"
    [forms_multicellular_complex]="yes no"
    [spore_formation]="yes no"
    [energy_source]="chemotroph phototroph"
    [electron_donor]="organotroph lithotroph"
    [carbon_source]="autotroph heterotroph"
    [biosafety_level]="1 2 3 4"
)

# BacDive columns copied to the selected genomes (as BACDIVE_<COLUMN>)
BACDIVE_ANNOTATION_TRAITS=(
    temperature_opt temperature_range ph_opt halophily_opt oxygen_tolerance cell_shape gram_stain
    spore_formation energy_source carbon_source
)
# Maximum number of distinct values per cell in the BACDIVE_* columns (0 = no limit)
MAX_LIST_ITEMS=20

# Prints the class that a user-typed value stands for (case-insensitive, common synonyms accepted).
# Returns 1 (prints nothing) if the value is not a class of the trait.
# Usage: canonical_class <trait> <value>
canonical_class() {
    local TRAIT="$1"
    local VALUE CLASS
    VALUE="$(printf '%s' "$2" | tr '[:upper:]' '[:lower:]' | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//; s/[[:space:]-]+/_/g')"

    case "${TRAIT}" in
        cell_shape)
            case "${VALUE}" in
                bacillus|bacilli|bacillary|rods|rod_shaped) VALUE="rod" ;;
                cocci|sphere|spherical)                     VALUE="coccus" ;;
                oval|ellipsoid|ellipsoidal|egg_shaped)      VALUE="ovoid" ;;
                vibrio|crescent|comma|bent)                 VALUE="curved" ;;
                helical|spirillum|spirilla)                 VALUE="spiral" ;;
                filamentous|hyphae)                         VALUE="filament" ;;
            esac
            ;;
        oxygen_tolerance)
            case "${VALUE}" in
                aerobic|aerobes|obligate_aerobe)                    VALUE="aerobe" ;;
                anaerobic|anaerobes|obligate_anaerobe)              VALUE="anaerobe" ;;
                microaerophilic|microaerobic|microaerophiles)       VALUE="microaerophile" ;;
                facultative_anaerobe|facultative_aerobe|facultatively_anaerobic) VALUE="facultative" ;;
            esac
            ;;
        gram_stain)
            case "${VALUE}" in
                pos|gram_positive|"+") VALUE="positive" ;;
                neg|gram_negative)     VALUE="negative" ;;
            esac
            ;;
        motility|pigmentation_production|forms_multicellular_complex|spore_formation)
            case "${VALUE}" in
                true|y)  VALUE="yes" ;;
                false|n) VALUE="no" ;;
            esac
            ;;
    esac

    for CLASS in ${TRAIT_CLASSES[${TRAIT}]}; do
        if [[ "${CLASS}" == "${VALUE}" ]]; then
            echo "${CLASS}"
            return 0
        fi
    done
    return 1
}

# Validates "<min>:<max>" (numbers; either side may be empty, not both; min <= max)
# Usage: validate_range <description> <range>
validate_range() {
    local LABEL="$1"
    local RANGE="$2"
    local MIN MAX VALUE
    local NUMBER_REGEX='^-?[0-9]+([.][0-9]+)?$'

    if [[ "${RANGE}" != *:* ]]; then
        die "Invalid range '${RANGE}' for ${LABEL}. Expected <min>:<max> (e.g., 45:80); either side may be left empty."
    fi
    MIN="${RANGE%%:*}"
    MAX="${RANGE#*:}"
    if [[ -z "${MIN}${MAX}" ]]; then
        die "Invalid range '${RANGE}' for ${LABEL}: give at least a minimum or a maximum."
    fi
    for VALUE in "${MIN}" "${MAX}"; do
        if [[ -n "${VALUE}" && ! "${VALUE}" =~ ${NUMBER_REGEX} ]]; then
            die "Invalid range '${RANGE}' for ${LABEL}: '${VALUE}' is not a number."
        fi
    done
    if [[ -n "${MIN}" && -n "${MAX}" ]] && awk -v a="${MIN}" -v b="${MAX}" 'BEGIN { exit !(a + 0 > b + 0) }'; then
        die "Invalid range '${RANGE}' for ${LABEL}: the minimum is greater than the maximum."
    fi
}

############ READING THE QUERY ################

Q_ACCESSION=()
Q_TAXID=()
Q_NCBI_TAXON=()
Q_GTDB_TAXON=()
Q_ISOLATION=()
Q_COUNTRY=()
Q_NCBI_RANGE=()      # "<column>=<min>:<max>"
Q_NCBI_MATCH=()      # "<column>=<text>"
Q_BACDIVE_RANGE=()
Q_BACDIVE_CLASS=()   # "<trait>=<class1>,<class2>" (canonical classes)
Q_GTDB_RANGE=()
Q_GTDB_MATCH=()
Q_ASSEMBLY_SOURCE="genbank"
Q_GENOME_TYPE="all"
declare -A SEEN_TRAIT=()

LINE_NUMBER=0
while IFS= read -r LINE || [[ -n "${LINE}" ]]; do
    LINE_NUMBER=$((LINE_NUMBER + 1))
    LINE="${LINE%$'\r'}"
    if [[ -z "${LINE//[[:space:]]/}" || "${LINE}" =~ ^[[:space:]]*# ]]; then
        continue
    fi
    if [[ "${LINE}" != *$'\t'* ]]; then
        die "Query line ${LINE_NUMBER} has no TAB between key and value: '${LINE}'"
    fi

    KEY="$(echo "${LINE%%$'\t'*}" | tr '[:lower:]' '[:upper:]' | tr -d '[:space:]')"
    VALUE="$(echo "${LINE#*$'\t'}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
    if [[ -z "${VALUE}" ]]; then
        die "Query line ${LINE_NUMBER} (${KEY}) has no value."
    fi

    case "${KEY}" in
        ACCESSION)
            VALUE="${VALUE^^}"
            if [[ ! "${VALUE}" =~ ^GC[AF]_[0-9]+([.][0-9]+)?$ ]]; then
                die "Invalid accession '${VALUE}' (expected e.g. GCA_000005845.2 or GCF_000005845)."
            fi
            Q_ACCESSION+=("${VALUE}")
            ;;
        TAXID)
            if [[ ! "${VALUE}" =~ ^[0-9]+$ ]]; then
                die "Invalid NCBI TaxID '${VALUE}' (expected a number)."
            fi
            Q_TAXID+=("${VALUE}")
            ;;
        NCBI_TAXON) Q_NCBI_TAXON+=("${VALUE}") ;;
        GTDB_TAXON) Q_GTDB_TAXON+=("${VALUE}") ;;
        ISOLATION_SOURCE) Q_ISOLATION+=("${VALUE}") ;;
        COUNTRY) Q_COUNTRY+=("${VALUE}") ;;
        ASSEMBLY_SOURCE)
            Q_ASSEMBLY_SOURCE="${VALUE,,}"
            if [[ "${Q_ASSEMBLY_SOURCE}" != "genbank" && "${Q_ASSEMBLY_SOURCE}" != "refseq" ]]; then
                die "Invalid assembly source '${VALUE}'. Must be 'genbank' or 'refseq'."
            fi
            ;;
        GENOME_TYPE)
            Q_GENOME_TYPE="${VALUE,,}"
            if [[ "${Q_GENOME_TYPE}" != "isolate_genome" && "${Q_GENOME_TYPE}" != "mag" && "${Q_GENOME_TYPE}" != "all" ]]; then
                die "Invalid genome type '${VALUE}'. Must be 'isolate_genome', 'mag', or 'all'."
            fi
            ;;
        NCBI_RANGE|BACDIVE_RANGE|GTDB_RANGE)
            if [[ "${VALUE}" != *=* || -z "${VALUE%%=*}" ]]; then
                die "Invalid ${KEY} '${VALUE}'. Expected <column>=<min>:<max>."
            fi
            validate_range "${KEY} ${VALUE%%=*}" "${VALUE#*=}"
            case "${KEY}" in
                NCBI_RANGE)    Q_NCBI_RANGE+=("${VALUE}") ;;
                BACDIVE_RANGE) Q_BACDIVE_RANGE+=("${VALUE}") ;;
                GTDB_RANGE)    Q_GTDB_RANGE+=("${VALUE}") ;;
            esac
            ;;
        NCBI_MATCH|GTDB_MATCH)
            if [[ "${VALUE}" != *=* || -z "${VALUE%%=*}" || -z "${VALUE#*=}" ]]; then
                die "Invalid ${KEY} '${VALUE}'. Expected <column>=<text>."
            fi
            if [[ "${KEY}" == "NCBI_MATCH" ]]; then Q_NCBI_MATCH+=("${VALUE}"); else Q_GTDB_MATCH+=("${VALUE}"); fi
            ;;
        BACDIVE_CLASS)
            if [[ "${VALUE}" != *=* ]]; then
                die "Invalid BACDIVE_CLASS '${VALUE}'. Expected <trait>=<class>[,<class>...] (e.g., gram_stain=negative)."
            fi
            TRAIT="$(printf '%s' "${VALUE%%=*}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
            if [[ ! "${TRAIT}" =~ ^[a-z0-9_]+$ || -z "${TRAIT_CLASSES[${TRAIT}]:-}" ]]; then
                die "Unknown trait '${VALUE%%=*}'. Traits divided in classes: ${TRAIT_ORDER[*]}."
            fi
            if [[ -n "${SEEN_TRAIT[${TRAIT}]:-}" ]]; then
                die "Trait '${TRAIT}' was given more than once. List several classes in one filter, separated by commas."
            fi
            SEEN_TRAIT["${TRAIT}"]=1
            CLASS_LIST=""
            IFS=',' read -ra ITEMS <<< "${VALUE#*=}"
            for ITEM in "${ITEMS[@]}"; do
                if [[ -z "${ITEM//[[:space:]]/}" ]]; then
                    continue
                fi
                if ! CANON="$(canonical_class "${TRAIT}" "${ITEM}")"; then
                    die "Invalid class '${ITEM}' for ${TRAIT}. Valid classes: ${TRAIT_CLASSES[${TRAIT}]// /, }."
                fi
                case ",${CLASS_LIST}," in
                    *",${CANON},"*) ;;
                    *) CLASS_LIST="${CLASS_LIST:+${CLASS_LIST},}${CANON}" ;;
                esac
            done
            if [[ -z "${CLASS_LIST}" ]]; then
                die "No class given for ${TRAIT}. Valid classes: ${TRAIT_CLASSES[${TRAIT}]// /, }."
            fi
            Q_BACDIVE_CLASS+=("${TRAIT}=${CLASS_LIST}")
            ;;
        *)
            die "Unknown key '${KEY}' in line ${LINE_NUMBER} of the query file."
            ;;
    esac
done < "${QUERY_FILE}"

HAS_BACDIVE_FILTER=$(( ${#Q_BACDIVE_RANGE[@]} + ${#Q_BACDIVE_CLASS[@]} > 0 ? 1 : 0 ))
HAS_GTDB_FILTER=$(( ${#Q_GTDB_TAXON[@]} + ${#Q_GTDB_RANGE[@]} + ${#Q_GTDB_MATCH[@]} > 0 ? 1 : 0 ))
N_NCBI_FILTERS=$(( ${#Q_ACCESSION[@]} + ${#Q_TAXID[@]} + ${#Q_NCBI_TAXON[@]} + ${#Q_ISOLATION[@]} + ${#Q_COUNTRY[@]} + ${#Q_NCBI_RANGE[@]} + ${#Q_NCBI_MATCH[@]} ))

# Without any filter, every genome of the NCBI tables would be selected (and downloaded)
if [[ "${HAS_BACDIVE_FILTER}" -eq 0 && "${HAS_GTDB_FILTER}" -eq 0 && "${N_NCBI_FILTERS}" -eq 0 && "${Q_GENOME_TYPE}" == "all" ]]; then
    die "The query has no selection filter (ASSEMBLY_SOURCE and 'GENOME_TYPE all' are not filters). Without filters every genome of the NCBI tables would be selected."
fi

if [[ "${VALIDATE_ONLY}" -eq 1 ]]; then
    exit 0
fi

############ REFERENCE TABLES ################

NCBI_DIR="${REF_DIR}/NCBI"
BACDIVE_TSV="$(pf_bacdive_table "${REF_DIR}/BACDIVE")"

NCBI_TABLES=()
for DOMAIN in "${PF_DOMAINS[@]}"; do
    TABLE="$(pf_ncbi_metadata_table "${NCBI_DIR}" "${DOMAIN}" "${Q_ASSEMBLY_SOURCE}")"
    if [[ ! -f "${TABLE}" ]]; then
        die "NCBI table not found: ${TABLE}. Was the NCBI reference database downloaded?"
    fi
    if [[ "${#NCBI_TABLES[@]}" -gt 0 && "$(head -n 1 "${TABLE}")" != "$(head -n 1 "${NCBI_TABLES[0]}")" ]]; then
        die "The columns of ${TABLE} differ from those of ${NCBI_TABLES[0]}."
    fi
    NCBI_TABLES+=("${TABLE}")
done
NCBI_HEADER_FILE="${NCBI_TABLES[0]}"

make_temp_dir WORK_DIR

# Common awk functions of the filters
AWK_LIB='
function is_empty(x)  { return x == "" || x == "NULL" || x == "NA" }
function is_number(x) { return x ~ /^[+-]?[0-9]+([.][0-9]+)?$/ }
# 1 if at least one of the values of the cell (separated by ";", ":" or "|") is inside [lo, hi]; empty bound = open
function any_in_range(x, lo, hi,   parts, n, i) {
    n = split(x, parts, /[;:|]/)
    for (i = 1; i <= n; i++) {
        if (!is_number(parts[i])) continue
        if (lo != "" && parts[i] + 0 < lo + 0) continue
        if (hi != "" && parts[i] + 0 > hi + 0) continue
        return 1
    }
    return 0
}
# Loads a spec file (<column index> <TAB> <arg1> <TAB> <arg2> <TAB> <label>) into the arrays s_col, s_a, s_b, s_label
function load_spec(file, s_col, s_a, s_b, s_label,   line, f, n) {
    if (file == "") return 0
    while ((getline line < file) > 0) {
        split(line, f, "\t")
        n++
        s_col[n] = f[1] + 0; s_a[n] = f[2]; s_b[n] = f[3]; s_label[n] = f[4]
    }
    close(file)
    return n + 0
}
# Loads a list file into a set (one entry per line)
function load_set(file, set,   line) {
    if (file == "") return
    while ((getline line < file) > 0) if (line != "") set[line] = 1
    close(file)
}
# GCA_000005845.2 -> 000005845.2 (GenBank and RefSeq assemblies of the same genome share the number)
function acc_number(acc,   a) { a = acc; sub(/^GC[AF]_/, "", a); return a }
function acc_base(acc,   a) { a = acc_number(acc); sub(/[.][0-9]+$/, "", a); return a }
'

# Writes a spec file for range/match filters: resolves each <column> name in the header of <table>.
# Usage: write_spec <kind: range|match> <table> <spec_file> <filter> [filter ...]
write_spec() {
    local KIND="$1"
    local TABLE="$2"
    local SPEC_FILE="$3"
    local FILTER COLUMN ARGUMENT INDEX
    shift 3

    : > "${SPEC_FILE}"
    for FILTER in "$@"; do
        COLUMN="${FILTER%%=*}"
        ARGUMENT="${FILTER#*=}"
        INDEX="$(get_column_index "${TABLE}" "${COLUMN}")"
        if [[ -z "${INDEX}" ]]; then
            die "Column '${COLUMN}' not found in ${TABLE}. Available columns: $(head -n 1 "${TABLE}" | tr '\t' ' ')"
        fi
        if [[ "${KIND}" == "range" ]]; then
            printf '%s\t%s\t%s\t%s\n' "${INDEX}" "${ARGUMENT%%:*}" "${ARGUMENT#*:}" "${COLUMN} ${ARGUMENT}" >> "${SPEC_FILE}"
        else
            printf '%s\t%s\t\t%s\n' "${INDEX}" "${ARGUMENT,,}" "${COLUMN}=${ARGUMENT}" >> "${SPEC_FILE}"
        fi
    done
}

# Logs the lines of a stats file written by awk ("<label><TAB><value>" per line)
log_stats() {
    local LABEL VALUE
    while IFS=$'\t' read -r LABEL VALUE; do
        log "message" "${LABEL}: ${VALUE}"
    done < "$1"
}

log "message" "Searching the reference database (NCBI ${Q_ASSEMBLY_SOURCE} tables)..."

############ STAGE 1: BACDIVE ################
#   Class filters: a strain passes if AT LEAST ONE of its classes (cells may hold several, separated by ";" or ":")
#   is among the classes requested for that trait. Range filters: a strain passes if ANY of its values is inside
#   the range. Strains with no value for a filtered trait are excluded. The NCBI TaxIDs of the strains that pass
#   are the link to NCBI.

BACDIVE_TAXIDS_FILE=""
if [[ "${HAS_BACDIVE_FILTER}" -eq 1 ]]; then
    if [[ ! -f "${BACDIVE_TSV}" ]]; then
        die "BacDive filters require the BacDive table, which was not found: ${BACDIVE_TSV}"
    fi
    log "message" "Stage 1: applying the phenotypic filters on BacDive..."

    BD_COL_TAXID="$(require_column "${BACDIVE_TSV}" ncbi_tax_id)"
    write_spec range "${BACDIVE_TSV}" "${WORK_DIR}/bacdive_range.spec" ${Q_BACDIVE_RANGE[@]+"${Q_BACDIVE_RANGE[@]}"}
    write_spec match "${BACDIVE_TSV}" "${WORK_DIR}/bacdive_class.spec" ${Q_BACDIVE_CLASS[@]+"${Q_BACDIVE_CLASS[@]}"}

    BACDIVE_TAXIDS_FILE="${WORK_DIR}/bacdive_taxids.txt"
    awk -F'\t' -v c_taxid="${BD_COL_TAXID}" \
        -v range_spec="${WORK_DIR}/bacdive_range.spec" -v class_spec="${WORK_DIR}/bacdive_class.spec" \
        -v stats="${WORK_DIR}/stage1.stats" "${AWK_LIB}"'
    BEGIN {
        nr = load_spec(range_spec, r_col, r_lo, r_hi, r_label)
        nc = load_spec(class_spec, c_col, c_classes, c_unused, c_label)
        for (i = 1; i <= nc; i++) {
            m = split(c_classes[i], want, ",")
            for (k = 1; k <= m; k++) class_want[i, want[k]] = 1
        }
    }
    function class_match(cell, i,   parts, n, j, tok) {
        if (is_empty(cell)) return 0
        n = split(tolower(cell), parts, /[;:]/)
        for (j = 1; j <= n; j++) {
            tok = parts[j]
            gsub(/^[ \t]+|[ \t]+$/, "", tok)
            if ((i, tok) in class_want) return 1
        }
        return 0
    }
    NR == 1 { next }
    {
        n_in++
        if (is_empty($c_taxid)) next
        n_linkable++
        for (i = 1; i <= nr; i++) { if (!any_in_range($(r_col[i]), r_lo[i], r_hi[i])) next; n_r[i]++ }
        for (i = 1; i <= nc; i++) { if (!class_match($(c_col[i]), i)) next; n_c[i]++ }
        # A strain can carry several TaxIDs joined with ";": each one is a candidate
        n_ids = split($c_taxid, ids, ";")
        for (j = 1; j <= n_ids; j++) if (!is_empty(ids[j])) print ids[j]
    }
    END {
        printf "  BacDive strains\t%d\n", n_in > stats
        printf "  ... with an NCBI TaxID (can be linked to NCBI)\t%d\n", n_linkable > stats
        for (i = 1; i <= nr; i++) printf "  ... after %s\t%d\n", r_label[i], n_r[i] > stats
        for (i = 1; i <= nc; i++) printf "  ... after %s\t%d\n", c_label[i], n_c[i] > stats
    }' "${BACDIVE_TSV}" | sort -u > "${BACDIVE_TAXIDS_FILE}"

    log_stats "${WORK_DIR}/stage1.stats"
    log "message" "Stage 1 complete: $(count_lines "${BACDIVE_TAXIDS_FILE}") distinct NCBI TaxIDs match the BacDive filters."
    if [[ ! -s "${BACDIVE_TAXIDS_FILE}" ]]; then
        log "warning" "No BacDive strain matched the phenotypic filters, so no genome can be selected."
    fi
else
    log "message" "Stage 1 skipped: no BacDive filter."
fi

############ STAGE 2: GTDB ################
#   GTDB_TAXON: the genome belongs to one of the taxa (any rank). GTDB_RANGE/GTDB_MATCH: GTDB metadata.
#   The result is a set of accession numbers (GCA_/GCF_ prefix and version removed), because GTDB lists each
#   genome under only one of its GenBank/RefSeq accessions.

GTDB_KEYS_FILE=""
if [[ "${HAS_GTDB_FILTER}" -eq 1 ]]; then
    log "message" "Stage 2: applying the GTDB filters..."
    GTDB_DIR="${REF_DIR}/GTDB"
    GTDB_TAX_TABLES=()
    GTDB_META_TABLES=()
    for DOMAIN in "${PF_DOMAINS[@]}"; do
        GTDB_TAX_TABLES+=("$(pf_gtdb_taxonomy_table "${GTDB_DIR}" "${DOMAIN}")")
        GTDB_META_TABLES+=("$(pf_gtdb_metadata_table "${GTDB_DIR}" "${DOMAIN}")")
    done

    GTDB_SETS=()

    if [[ "${#Q_GTDB_TAXON[@]}" -gt 0 ]]; then
        for TABLE in "${GTDB_TAX_TABLES[@]}"; do
            [[ -f "${TABLE}" ]] || die "GTDB taxonomy table not found: ${TABLE}. Was the GTDB reference database downloaded?"
        done
        printf '%s\n' "${Q_GTDB_TAXON[@]}" | awk '{ v = tolower($0); sub(/^[a-z]__/, "", v); print v }' > "${WORK_DIR}/gtdb_taxa.txt"
        awk -F'\t' -v taxa_file="${WORK_DIR}/gtdb_taxa.txt" "${AWK_LIB}"'
            BEGIN { load_set(taxa_file, want) }
            FNR == 1 { for (i = 1; i <= NF; i++) if ($i ~ /^GTDB_(DOMAIN|PHYLUM|CLASS|ORDER|FAMILY|GENUS|SPECIES)$/) rank_col[i] = 1; next }
            { for (i in rank_col) if (tolower($i) in want) { print acc_base($1); next } }
        ' "${GTDB_TAX_TABLES[@]}" | sort -u > "${WORK_DIR}/gtdb_taxon_keys.txt"
        log "message" "  GTDB genomes in the requested taxa (${Q_GTDB_TAXON[*]}): $(count_lines "${WORK_DIR}/gtdb_taxon_keys.txt")"
        GTDB_SETS+=("${WORK_DIR}/gtdb_taxon_keys.txt")
    fi

    if [[ $(( ${#Q_GTDB_RANGE[@]} + ${#Q_GTDB_MATCH[@]} )) -gt 0 ]]; then
        for TABLE in "${GTDB_META_TABLES[@]}"; do
            [[ -f "${TABLE}" ]] || die "GTDB metadata table not found: ${TABLE}. Was the GTDB reference database downloaded?"
        done
        write_spec range "${GTDB_META_TABLES[0]}" "${WORK_DIR}/gtdb_range.spec" ${Q_GTDB_RANGE[@]+"${Q_GTDB_RANGE[@]}"}
        write_spec match "${GTDB_META_TABLES[0]}" "${WORK_DIR}/gtdb_match.spec" ${Q_GTDB_MATCH[@]+"${Q_GTDB_MATCH[@]}"}
        awk -F'\t' -v range_spec="${WORK_DIR}/gtdb_range.spec" -v match_spec="${WORK_DIR}/gtdb_match.spec" \
            -v stats="${WORK_DIR}/stage2.stats" "${AWK_LIB}"'
            BEGIN {
                nr = load_spec(range_spec, r_col, r_lo, r_hi, r_label)
                nm = load_spec(match_spec, m_col, m_text, m_unused, m_label)
            }
            FNR == 1 { next }
            {
                n_in++
                for (i = 1; i <= nr; i++) { if (!any_in_range($(r_col[i]), r_lo[i], r_hi[i])) next; n_r[i]++ }
                for (i = 1; i <= nm; i++) { if (is_empty($(m_col[i])) || index(tolower($(m_col[i])), m_text[i]) == 0) next; n_m[i]++ }
                print acc_base($1)
            }
            END {
                printf "  GTDB genomes (metadata)\t%d\n", n_in > stats
                for (i = 1; i <= nr; i++) printf "  ... after %s\t%d\n", r_label[i], n_r[i] > stats
                for (i = 1; i <= nm; i++) printf "  ... after %s\t%d\n", m_label[i], n_m[i] > stats
            }
        ' "${GTDB_META_TABLES[@]}" | sort -u > "${WORK_DIR}/gtdb_metadata_keys.txt"
        log_stats "${WORK_DIR}/stage2.stats"
        GTDB_SETS+=("${WORK_DIR}/gtdb_metadata_keys.txt")
    fi

    # Intersection of the active GTDB filters
    GTDB_KEYS_FILE="${WORK_DIR}/gtdb_keys.txt"
    cp "${GTDB_SETS[0]}" "${GTDB_KEYS_FILE}"
    for SET in "${GTDB_SETS[@]:1}"; do
        comm -12 "${GTDB_KEYS_FILE}" "${SET}" > "${GTDB_KEYS_FILE}.tmp"
        mv "${GTDB_KEYS_FILE}.tmp" "${GTDB_KEYS_FILE}"
    done
    log "message" "Stage 2 complete: $(count_lines "${GTDB_KEYS_FILE}") GTDB genomes match the GTDB filters."
else
    log "message" "Stage 2 skipped: no GTDB filter."
fi

############ STAGE 3: NCBI ################

log "message" "Stage 3: querying the NCBI metadata tables..."

C_ACC="$(require_column "${NCBI_HEADER_FILE}" ASSEMBLY_ACCESSION)"
C_TAXID="$(require_column "${NCBI_HEADER_FILE}" ORGANISM_TAXONOMIC_ID)"
C_ORG="$(require_column "${NCBI_HEADER_FILE}" ORGANISM_NAME)"
C_ISO="$(require_column "${NCBI_HEADER_FILE}" ASSEMBLY_BIOSAMPLE_ISOLATION_SOURCE)"
C_COUNTRY="$(require_column "${NCBI_HEADER_FILE}" ASSEMBLY_BIOSAMPLE_GEOGRAPHIC_LOCATION)"
C_GTYPE="$(require_column "${NCBI_HEADER_FILE}" GENOME_TYPE)"
C_SPECIES_TAXID="$(require_column "${NCBI_HEADER_FILE}" NCBI_SPECIES_TAXID)"
C_STRAIN_TAXID="$(require_column "${NCBI_HEADER_FILE}" NCBI_STRAIN_TAXID)"

# Columns of the NCBI lineage (names and TaxIDs of each rank)
LINEAGE_NAME_COLS="${C_ORG}"
LINEAGE_TAXID_COLS="${C_TAXID}"
for RANK in DOMAIN PHYLUM CLASS ORDER FAMILY GENUS SPECIES STRAIN; do
    LINEAGE_NAME_COLS="${LINEAGE_NAME_COLS},$(require_column "${NCBI_HEADER_FILE}" "NCBI_${RANK}")"
    LINEAGE_TAXID_COLS="${LINEAGE_TAXID_COLS},$(require_column "${NCBI_HEADER_FILE}" "NCBI_${RANK}_TAXID")"
done

# Lists of the OR filters (empty path = filter off)
write_list() {
    if [[ "$#" -gt 1 ]]; then
        printf '%s\n' "${@:2}" > "$1"
        echo "$1"
    fi
}
ACC_LIST="$(write_list "${WORK_DIR}/accessions.txt" ${Q_ACCESSION[@]+"${Q_ACCESSION[@]}"})"
TAXID_LIST="$(write_list "${WORK_DIR}/taxids.txt" ${Q_TAXID[@]+"${Q_TAXID[@]}"})"
TAXON_LIST="$(write_list "${WORK_DIR}/ncbi_taxa.txt" ${Q_NCBI_TAXON[@]+"${Q_NCBI_TAXON[@]}"})"
ISO_LIST="$(write_list "${WORK_DIR}/isolation.txt" ${Q_ISOLATION[@]+"${Q_ISOLATION[@]}"})"
COUNTRY_LIST="$(write_list "${WORK_DIR}/countries.txt" ${Q_COUNTRY[@]+"${Q_COUNTRY[@]}"})"
write_spec range "${NCBI_HEADER_FILE}" "${WORK_DIR}/ncbi_range.spec" ${Q_NCBI_RANGE[@]+"${Q_NCBI_RANGE[@]}"}
write_spec match "${NCBI_HEADER_FILE}" "${WORK_DIR}/ncbi_match.spec" ${Q_NCBI_MATCH[@]+"${Q_NCBI_MATCH[@]}"}

SELECTED_TMP="${WORK_DIR}/selected.tsv"

awk -F'\t' -v OFS='\t' \
    -v c_acc="${C_ACC}" -v c_taxid="${C_TAXID}" -v c_iso="${C_ISO}" -v c_country="${C_COUNTRY}" \
    -v c_gtype="${C_GTYPE}" -v c_species_taxid="${C_SPECIES_TAXID}" -v c_strain_taxid="${C_STRAIN_TAXID}" \
    -v name_cols="${LINEAGE_NAME_COLS}" -v taxid_cols="${LINEAGE_TAXID_COLS}" \
    -v bacdive_file="${BACDIVE_TAXIDS_FILE}" -v gtdb_file="${GTDB_KEYS_FILE}" \
    -v acc_file="${ACC_LIST}" -v taxid_file="${TAXID_LIST}" -v taxon_file="${TAXON_LIST}" \
    -v iso_file="${ISO_LIST}" -v country_file="${COUNTRY_LIST}" \
    -v range_spec="${WORK_DIR}/ncbi_range.spec" -v match_spec="${WORK_DIR}/ncbi_match.spec" \
    -v gtype="${Q_GENOME_TYPE}" -v stats="${WORK_DIR}/stage3.stats" "${AWK_LIB}"'
function norm(s) { s = tolower(s); gsub(/[^a-z0-9]+/, " ", s); sub(/ +$/, "", s); return " " s }
function load_list(file, list, as_words,   line, n) {
    if (file == "") return 0
    while ((getline line < file) > 0) if (line != "") list[++n] = (as_words ? norm(line) : tolower(line))
    close(file)
    return n + 0
}
function contains_any(text, list, n,   i) {
    for (i = 1; i <= n; i++) if (index(text, list[i]) > 0) return 1
    return 0
}
BEGIN {
    has_bd = (bacdive_file != ""); has_gtdb = (gtdb_file != "")
    has_acc = (acc_file != ""); has_taxid = (taxid_file != ""); has_taxon = (taxon_file != "")
    has_iso = (iso_file != ""); has_country = (country_file != "")
    has_gtype = (gtype != "all")

    load_set(bacdive_file, bd_taxids)
    load_set(gtdb_file, gtdb_keys)
    load_set(taxid_file, wanted_taxids)
    if (has_acc) {
        while ((getline a < acc_file) > 0) wanted_acc[acc_number(a)] = 1
        close(acc_file)
    }
    if (has_taxon) { while ((getline a < taxon_file) > 0) wanted_taxon[tolower(a)] = 1; close(taxon_file) }
    n_iso = load_list(iso_file, iso_list, 0)
    n_country = load_list(country_file, country_list, 1)
    nr = load_spec(range_spec, r_col, r_lo, r_hi, r_label)
    nm = load_spec(match_spec, m_col, m_text, m_unused, m_label)
    n_name_cols = split(name_cols, name_col, ",")
    n_taxid_cols = split(taxid_cols, taxid_col, ",")
}
FNR == 1 { if (NR == 1) print; next }
{
    n_in++

    # 1. BacDive TaxIDs (stage 1): organism, strain or species TaxID of the assembly
    if (has_bd && !(($c_taxid in bd_taxids) || ($c_strain_taxid in bd_taxids) || ($c_species_taxid in bd_taxids))) next
    n_bd++

    # 2. GTDB genomes (stage 2)
    if (has_gtdb && !(acc_base($c_acc) in gtdb_keys)) next
    n_gtdb++

    # 3. Accessions (version optional; GenBank/RefSeq pairs match each other)
    if (has_acc && !((acc_number($c_acc) in wanted_acc) || (acc_base($c_acc) in wanted_acc))) next
    n_acc++

    # 4. TaxID of the organism or of any rank of its lineage
    if (has_taxid) {
        ok = 0
        for (i = 1; i <= n_taxid_cols; i++) if ($(taxid_col[i]) in wanted_taxids) { ok = 1; break }
        if (!ok) next
    }
    n_taxid++

    # 5. NCBI taxon names
    if (has_taxon) {
        ok = 0
        for (i = 1; i <= n_name_cols; i++) if (tolower($(name_col[i])) in wanted_taxon) { ok = 1; break }
        if (!ok) next
    }
    n_taxon++

    # 6. Isolation source (substring) and 7. country (start of a word)
    if (has_iso && (is_empty($c_iso) || !contains_any(tolower($c_iso), iso_list, n_iso))) next
    n_isolation++
    if (has_country && (is_empty($c_country) || !contains_any(norm($c_country), country_list, n_country))) next
    n_country_ok++

    # 8. Genome type
    if (has_gtype && tolower($c_gtype) != gtype) next
    n_gtype++

    # 9. Ranges and text matches on any NCBI column
    for (i = 1; i <= nr; i++) { if (!any_in_range($(r_col[i]), r_lo[i], r_hi[i])) next; n_r[i]++ }
    for (i = 1; i <= nm; i++) { if (is_empty($(m_col[i])) || index(tolower($(m_col[i])), m_text[i]) == 0) next; n_m[i]++ }

    print
}
END {
    printf "  NCBI genome assemblies searched\t%d\n", n_in > stats
    if (has_bd)      printf "  ... linked to the BacDive TaxIDs (stage 1)\t%d\n", n_bd > stats
    if (has_gtdb)    printf "  ... in the GTDB genomes (stage 2)\t%d\n", n_gtdb > stats
    if (has_acc)     printf "  ... after accession filter\t%d\n", n_acc > stats
    if (has_taxid)   printf "  ... after TaxID filter\t%d\n", n_taxid > stats
    if (has_taxon)   printf "  ... after NCBI taxon filter\t%d\n", n_taxon > stats
    if (has_iso)     printf "  ... after isolation source filter\t%d\n", n_isolation > stats
    if (has_country) printf "  ... after country filter\t%d\n", n_country_ok > stats
    if (has_gtype)   printf "  ... after genome type filter\t%d\n", n_gtype > stats
    for (i = 1; i <= nr; i++) printf "  ... after %s\t%d\n", r_label[i], n_r[i] > stats
    for (i = 1; i <= nm; i++) printf "  ... after %s\t%d\n", m_label[i], n_m[i] > stats
}' "${NCBI_TABLES[@]}" > "${SELECTED_TMP}"

log_stats "${WORK_DIR}/stage3.stats"
TOTAL_SELECTED="$(count_records "${SELECTED_TMP}")"
log "message" "Stage 3 complete: ${TOTAL_SELECTED} genome assemblies selected."

############ BACDIVE ANNOTATION ################
#   The BacDive traits of each selected assembly are added as BACDIVE_* columns. An assembly is linked to the
#   BacDive strains of its organism TaxID (BACDIVE_LINK_LEVEL = organism) or, if there are none, to those of its
#   species TaxID (BACDIVE_LINK_LEVEL = species). Distinct values of the linked strains are joined with "|".

if [[ -f "${BACDIVE_TSV}" && "${TOTAL_SELECTED}" -gt 0 ]]; then
    log "message" "Adding the BacDive traits of the selected genomes..."
    BD_SPEC=""
    for TRAIT in bacdive_id ${BACDIVE_ANNOTATION_TRAITS[@]+"${BACDIVE_ANNOTATION_TRAITS[@]}"}; do
        INDEX="$(get_column_index "${BACDIVE_TSV}" "${TRAIT}")"
        if [[ -n "${INDEX}" ]]; then
            BD_SPEC="${BD_SPEC:+${BD_SPEC},}${INDEX}=${TRAIT}"
        fi
    done

    awk -F'\t' -v OFS='\t' -v bd_file="${BACDIVE_TSV}" -v bd_taxid_col="$(require_column "${BACDIVE_TSV}" ncbi_tax_id)" \
        -v spec="${BD_SPEC}" -v c_taxid="${C_TAXID}" -v c_species_taxid="${C_SPECIES_TAXID}" \
        -v max_items="${MAX_LIST_ITEMS}" '
    function is_empty(x) { return x == "" || x == "NULL" || x == "NA" }
    function add(t, k, v,   key) {
        if (is_empty(v)) return
        key = t SUBSEP k SUBSEP v
        if (key in seen) return
        seen[key] = 1
        n = ++cnt[t, k]
        if (max_items + 0 == 0 || n <= max_items + 0) {
            if ((t, k) in val) val[t, k] = val[t, k] "|" v
            else val[t, k] = v
        }
    }
    function get(t, k,   s) {
        if (!((t, k) in val)) return "NULL"
        s = val[t, k]
        if (max_items + 0 > 0 && cnt[t, k] > max_items + 0) s = s "|...(+" (cnt[t, k] - max_items) " more)"
        return s
    }
    BEGIN {
        ns = split(spec, sp, ",")
        for (k = 1; k <= ns; k++) { split(sp[k], p, "="); col[k] = p[1] + 0; name[k] = toupper(p[2]) }
        while ((getline line < bd_file) > 0) {
            if (++rows == 1) continue
            split(line, f, "\t")
            n_ids = split(f[bd_taxid_col], ids, ";")
            for (j = 1; j <= n_ids; j++) {
                if (is_empty(ids[j])) continue
                linked[ids[j]] = 1
                for (k = 1; k <= ns; k++) add(ids[j], k, f[col[k]])
            }
        }
        close(bd_file)
        out = "ASSEMBLY_ACCESSION" OFS "BACDIVE_LINK_LEVEL"
        for (k = 1; k <= ns; k++) out = out OFS (name[k] == "BACDIVE_ID" ? "BACDIVE_IDS" : "BACDIVE_" name[k])
        print out
    }
    NR == 1 { next }
    {
        if ($c_taxid in linked)              { t = $c_taxid; level = "organism" }
        else if ($c_species_taxid in linked) { t = $c_species_taxid; level = "species" }
        else                                 { t = ""; level = "NULL" }
        out = $1 OFS level
        for (k = 1; k <= ns; k++) out = out OFS (t == "" ? "NULL" : get(t, k))
        print out
    }' "${SELECTED_TMP}" > "${WORK_DIR}/bacdive_annotation.tsv"

    add_metadata_columns "${SELECTED_TMP}" "${WORK_DIR}/bacdive_annotation.tsv"
    log "message" "  Selected genomes linked to BacDive: $(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "BACDIVE_LINK_LEVEL") c = i; next } $c != "NULL"' "${SELECTED_TMP}" | wc -l)"
fi

mkdir -p "$(dirname "${OUTPUT_TSV}")"
mv "${SELECTED_TMP}" "${OUTPUT_TSV}"
log "message" "Selected genome metadata (one row per assembly) saved to: ${OUTPUT_TSV}"
