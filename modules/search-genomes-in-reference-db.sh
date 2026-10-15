#!/usr/bin/env bash
#
# search-genomes-in-reference-db.sh - Prokfinder module 1: genome search in the reference database
#
# Reads a query file (genome-selection-query.txt) and selects, in the NCBI tables, the genome assemblies that
# satisfy every filter. Filters can use the metadata of NCBI, the phenotypic traits of BacDive and the
# taxonomy/metadata of GTDB:
#
#   Stage 1 (BacDive) -> NCBI TaxIDs whose BacDive traits satisfy each BACDIVE_* filter. The BacDive table used has
#                        one row per NCBI TaxID (the strains of the same TaxID merged, see
#                        download-reference-db-bacdive.sh). A TaxID of rank species is linked to the assemblies of
#                        that species (NCBI_SPECIES_TAXID); a TaxID below the species rank (strain, subspecies...) to
#                        the assemblies of that organism or strain (exact TaxID); TaxIDs above the species rank
#                        (genus, family...) are not linked. An assembly linked to a TaxID below the species rank
#                        uses its row; otherwise the row of its species.
#   Stage 2 (GTDB)    -> accessions whose GTDB taxonomy/metadata satisfy each GTDB_* filter
#   Stage 3 (NCBI)    -> each NCBI assembly is checked against every filter; the assemblies that satisfy all of
#                        them are selected
#   Annotation        -> BacDive traits of the selected assemblies are added as BACDIVE_* columns
#
# QUERY FILE: one filter per line, "<KEY><TAB><VALUE>". Empty lines and lines starting with "#" are ignored.
# Lines with the same key in the first group are combined with OR; everything else is combined with AND.
#
#   OR within the key:
#     ACCESSION         GCA_/GCF_ accession (version optional). GCA_x.1 also matches its paired GCF_x.1
#     TAXID             NCBI TaxID of the organism or of any rank of its lineage (descendant taxa included)
#     NCBI_TAXON        name of a taxon of the NCBI lineage (any rank, exact match, case-insensitive)
#     TAXON_NAME        name of an NCBI taxon, converted to its TaxID with "taxonkit name2taxid" (Bacteria and
#                       Archaea only); then the same as TAXID. A name without TaxID is an error
#     GTDB_TAXON        name of a taxon of the GTDB lineage (any rank; "g__Thermus" or "Thermus")
#     ISOLATION_SOURCE  text contained in the isolation source (case-insensitive)
#     COUNTRY           country/place at the start of a word of the geographic location (case-insensitive)
#   Single value:
#     ASSEMBLY_SOURCE   genbank (default) or refseq: NCBI table that is searched
#     GENOME_TYPE       isolate_genome, mag or all (default)
#   AND (each line is one more condition):
#     NCBI_RANGE        <column>=<min>:<max>    numeric column of the NCBI table (either side may be empty)
#     NCBI_MATCH        <column>=<text>         text contained in a column of the NCBI table
#     BACDIVE_RANGE     <column>=<min>:<max>    numeric trait of BacDive (the value of the TaxID inside the range)
#     BACDIVE_CLASS     <trait>=<class>[,<class>...]   trait divided in classes (any of the classes was recorded)
#     GTDB_RANGE        <column>=<min>:<max>    numeric column of the GTDB metadata
#     GTDB_MATCH        <column>=<text>         text contained in a column of the GTDB metadata
#
# Usage:
#   search-genomes-in-reference-db.sh -q <query_file> --validate-only [--taxdump-dir <dir>]
#   search-genomes-in-reference-db.sh -q <query_file> -r <reference_db_dir> -o <selected_metadata.tsv> [-f <filters.tsv>]
#
# Output:
#   <selected_metadata.tsv>  one row per selected assembly: NCBI columns, the *_MATCH columns (value that matched
#                            each filter with several values, e.g. ISOLATION_SOURCE_KEYWORD_MATCH) and the
#                            BACDIVE_* columns
#   <filters.tsv>            every assembly that satisfies AT LEAST ONE filter: SELECTED (TRUE when it satisfies all
#                            of them), one SELECT_BY_<FILTER> column (TRUE/FALSE) per filter and the *_MATCH columns
#                            (default: genome-selection-filters.tsv next to <selected_metadata.tsv>)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -q <query_file> [-r <reference_db_dir> -o <output_tsv> [-f <filters_tsv>]] [--validate-only]

  -q <file>             Query file (see the header of this script for the format)
  -r <dir>              Reference databases folder (~/prokfinder-refdb)
  -o <file>             Output table with the selected genome assemblies
  -f <file>             Table of the assemblies that satisfy at least one filter, with TRUE/FALSE per filter
                        (default: genome-selection-filters.tsv in the folder of -o)
  --validate-only       Only check the query file (no reference table is read)
  --taxdump-dir <dir>   NCBI taxdump used to check the TAXON_NAME names with --validate-only (when it is not
                        given, the names are checked only if taxonkit has its default taxdump, ~/.taxonkit)
  -h                    Show this help message
EOF
    exit "${1:-1}"
}

QUERY_FILE=""
REF_DIR=""
OUTPUT_TSV=""
FILTERS_TSV=""
TAXDUMP_DIR=""
VALIDATE_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -q) QUERY_FILE="$2"; shift 2 ;;
        -r) REF_DIR="${2%/}"; shift 2 ;;
        -o) OUTPUT_TSV="$2"; shift 2 ;;
        -f) FILTERS_TSV="$2"; shift 2 ;;
        --validate-only) VALIDATE_ONLY=1; shift ;;
        --taxdump-dir) TAXDUMP_DIR="${2%/}"; shift 2 ;;
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
FILTERS_TSV="${FILTERS_TSV:-$(dirname "${OUTPUT_TSV:-.}")/genome-selection-filters.tsv}"

if [[ "${VALIDATE_ONLY}" -eq 0 ]]; then
    setup_module_log "$(dirname "${OUTPUT_TSV}")/search-genomes-in-reference-db.log"
fi

make_temp_dir WORK_DIR

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
Q_TAXON_NAME=()
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
BAD_ACCESSIONS=()

LINE_NUMBER=0
while IFS= read -r LINE || [[ -n "${LINE}" ]]; do
    LINE_NUMBER=$((LINE_NUMBER + 1))
    LINE="${LINE//$'\r'/}"
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
            if [[ "${VALUE}" =~ ^GC[AF]_[0-9]+([.][0-9]+)?$ ]]; then
                Q_ACCESSION+=("${VALUE}")
            else
                BAD_ACCESSIONS+=("${VALUE}")
            fi
            ;;
        TAXID)
            if [[ ! "${VALUE}" =~ ^[0-9]+$ ]]; then
                die "Invalid NCBI TaxID '${VALUE}' (expected a number)."
            fi
            Q_TAXID+=("${VALUE}")
            ;;
        NCBI_TAXON) Q_NCBI_TAXON+=("${VALUE}") ;;
        TAXON_NAME) Q_TAXON_NAME+=("${VALUE}") ;;
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

# Every invalid accession is reported at once
if [[ "${#BAD_ACCESSIONS[@]}" -gt 0 ]]; then
    LIST="$(printf '%s\n' "${BAD_ACCESSIONS[@]}" | head -n 20 | paste -sd ' ')"
    MORE=$(( ${#BAD_ACCESSIONS[@]} > 20 ? ${#BAD_ACCESSIONS[@]} - 20 : 0 ))
    die "${#BAD_ACCESSIONS[@]} invalid accession(s): ${LIST}$( (( MORE > 0 )) && echo " (and ${MORE} more)"). An accession is GCA_ or GCF_ followed by the number and, optionally, the version (e.g., GCA_000005845.2 or GCF_000005845)."
fi

HAS_BACDIVE_FILTER=$(( ${#Q_BACDIVE_RANGE[@]} + ${#Q_BACDIVE_CLASS[@]} > 0 ? 1 : 0 ))
HAS_GTDB_FILTER=$(( ${#Q_GTDB_TAXON[@]} + ${#Q_GTDB_RANGE[@]} + ${#Q_GTDB_MATCH[@]} > 0 ? 1 : 0 ))
N_NCBI_FILTERS=$(( ${#Q_ACCESSION[@]} + ${#Q_TAXID[@]} + ${#Q_NCBI_TAXON[@]} + ${#Q_TAXON_NAME[@]} + ${#Q_ISOLATION[@]} + ${#Q_COUNTRY[@]} + ${#Q_NCBI_RANGE[@]} + ${#Q_NCBI_MATCH[@]} ))

# Without any filter, every genome of the NCBI tables would be selected (and downloaded)
if [[ "${HAS_BACDIVE_FILTER}" -eq 0 && "${HAS_GTDB_FILTER}" -eq 0 && "${N_NCBI_FILTERS}" -eq 0 && "${Q_GENOME_TYPE}" == "all" ]]; then
    die "The query has no selection filter (ASSEMBLY_SOURCE and 'GENOME_TYPE all' are not filters). Without filters every genome of the NCBI tables would be selected."
fi

############ TAXON NAMES -> NCBI TAXIDS ################

# Converts the TAXON_NAME names into NCBI TaxIDs of Bacteria or Archaea (taxonkit name2taxid). Stops with the list of
# the names that have no such TaxID. A name with several TaxIDs (homonyms) keeps all of them.
# Usage: resolve_taxon_names <taxdump_dir, or "" for the taxonkit default> <output: taxid TAB name>
resolve_taxon_names() {
    local TAXDUMP="$1"
    local OUT="$2"
    local -a ARGS=()
    local MISSING
    if [[ -n "${TAXDUMP}" ]]; then
        ARGS=(--data-dir "${TAXDUMP}")
    fi

    # name2taxid prints "<name> <TAB> <taxid>" (one line per TaxID; empty TaxID when the name is unknown)
    printf '%s\n' "${Q_TAXON_NAME[@]}" \
        | taxonkit name2taxid ${ARGS[@]+"${ARGS[@]}"} > "${WORK_DIR}/name2taxid.tsv" 2> "${WORK_DIR}/name2taxid.err" \
        || die "taxonkit name2taxid failed: $(head -n 3 "${WORK_DIR}/name2taxid.err" | tr '\n' ' ')"
    awk -F'\t' '$2 ~ /^[0-9]+$/ { print $2 }' "${WORK_DIR}/name2taxid.tsv" | sort -u \
        | taxonkit lineage -t ${ARGS[@]+"${ARGS[@]}"} > "${WORK_DIR}/name2taxid.lineage" 2>> "${WORK_DIR}/name2taxid.err" \
        || die "taxonkit lineage failed: $(head -n 3 "${WORK_DIR}/name2taxid.err" | tr '\n' ' ')"

    # Only TaxIDs whose lineage contains Bacteria (2) or Archaea (2157)
    awk -F'\t' -v OFS='\t' -v lineage="${WORK_DIR}/name2taxid.lineage" '
        BEGIN { while ((getline line < lineage) > 0) { split(line, f, "\t"); if ((";" f[3] ";") ~ /;(2|2157);/) prok[f[1]] = 1 } }
        $2 ~ /^[0-9]+$/ && ($2 in prok) { print $2, $1 }
    ' "${WORK_DIR}/name2taxid.tsv" | sort -u > "${OUT}"

    MISSING="$(printf '%s\n' "${Q_TAXON_NAME[@]}" | awk -F'\t' -v found="${OUT}" '
        BEGIN { while ((getline line < found) > 0) { split(line, f, "\t"); ok[tolower(f[2])] = 1 } }
        !(tolower($0) in ok) { printf "%s\"%s\"", (n++ ? ", " : ""), $0 }')"
    if [[ -n "${MISSING}" ]]; then
        die "Taxon name(s) not found in the NCBI taxonomy of Bacteria and Archaea: ${MISSING}. Check the spelling (scientific names, e.g., Thermus or Thermus aquaticus)."
    fi
}

if [[ "${VALIDATE_ONLY}" -eq 1 ]]; then
    # The taxon names are checked now only if taxonkit and a taxdump are available here; otherwise, during the search
    if [[ "${#Q_TAXON_NAME[@]}" -gt 0 ]] && command -v taxonkit >/dev/null 2>&1; then
        if [[ -n "${TAXDUMP_DIR}" && -f "${TAXDUMP_DIR}/names.dmp" ]]; then
            resolve_taxon_names "${TAXDUMP_DIR}" "${WORK_DIR}/taxon_names.tsv"
        elif [[ -f "${HOME}/.taxonkit/names.dmp" ]]; then
            resolve_taxon_names "" "${WORK_DIR}/taxon_names.tsv"
        fi
    fi
    exit 0
fi

############ REFERENCE TABLES ################

NCBI_DIR="${REF_DIR}/NCBI"
BACDIVE_TSV="$(pf_bacdive_species_table "${REF_DIR}/BACDIVE")"

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
# Ranks above species: BacDive TaxIDs of these ranks are not linked to genomes
function above_species(r) {
    return r ~ /^(superkingdom|domain|realm|kingdom|subkingdom|phylum|subphylum|superclass|class|subclass|infraclass|superorder|order|suborder|infraorder|parvorder|superfamily|family|subfamily|tribe|subtribe|genus|subgenus|section|subsection|series|species group|species subgroup)$/
}
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

# Every filter becomes a SELECT_BY_<NAME> column of the filters table. The filters are registered in the order in
# which stage 3 checks them (the names are made unique with a number when a filter is repeated).
FILTER_NAMES_FILE="${WORK_DIR}/filter_names.tsv"
: > "${FILTER_NAMES_FILE}"
declare -A USED_FILTER_NAME=()
# Usage: add_filter_name <name> <label>
add_filter_name() {
    local NAME N=2
    NAME="SELECT_BY_$(printf '%s' "$1" | tr '[:lower:]' '[:upper:]' | sed -E 's/[^A-Z0-9]+/_/g; s/^_+|_+$//g')"
    if [[ -n "${USED_FILTER_NAME[${NAME}]:-}" ]]; then
        while [[ -n "${USED_FILTER_NAME[${NAME}_${N}]:-}" ]]; do N=$((N + 1)); done
        NAME="${NAME}_${N}"
    fi
    USED_FILTER_NAME["${NAME}"]=1
    printf '%s\t%s\n' "${NAME}" "$2" >> "${FILTER_NAMES_FILE}"
}

# Name of the SELECT_BY column of a BacDive filter (named after the option of main-prokfinder.sh when there is one;
# otherwise the trait for a class filter, BACDIVE_<column> for a range filter)
# Usage: bacdive_filter_name <range|class> <column>
bacdive_filter_name() {
    if [[ "$1" == "class" && "${2,,}" != "cell_shape" && "${2,,}" != "oxygen_tolerance" ]]; then
        echo "$2"
        return 0
    fi
    shift
    case "${1,,}" in
        temperature_opt)        echo "OPT_TEMP" ;;
        ph_opt)                 echo "OPT_PH" ;;
        halophily_opt)          echo "OPT_SALT" ;;
        temperature_preference) echo "TEMP_GROWTH" ;;
        ph_preference)          echo "PH_GROWTH" ;;
        halophily_preference)   echo "SALT_GROWTH" ;;
        cell_shape)             echo "CELL_SHAPE" ;;
        oxygen_tolerance)       echo "OXYGEN" ;;
        *)                      echo "BACDIVE_$1" ;;
    esac
}

log "message" "Searching the reference database (NCBI ${Q_ASSEMBLY_SOURCE} tables)..."

############ STAGE 1: BACDIVE ################
#   One row per NCBI TaxID (the BacDive strains of the same TaxID merged). Range filters: the value of the TaxID
#   (mean of its strains) is inside the range. Class filters: at least one of the requested classes was recorded for
#   the TaxID. TaxIDs with no value for a filtered trait do not pass. Output: one line per TaxID that can be linked to
#   genomes: <taxid> <TAB> <species|exact> <TAB> <rank> <TAB> <one 0/1 flag per BacDive filter>.

BD_LINKS_FILE=""
if [[ "${HAS_BACDIVE_FILTER}" -eq 1 ]]; then
    if [[ ! -f "${BACDIVE_TSV}" ]]; then
        die "BacDive filters require the BacDive TaxID table, which was not found: ${BACDIVE_TSV}. Run Prokfinder again (the BacDive tables are rebuilt automatically from the files already downloaded)."
    fi
    log "message" "Stage 1: applying the phenotypic filters on BacDive (one row per NCBI TaxID, strains merged)..."

    BD_COL_TAXID="$(require_column "${BACDIVE_TSV}" ncbi_tax_id)"
    BD_COL_RANK="$(require_column "${BACDIVE_TSV}" ncbi_tax_id_rank)"
    write_spec range "${BACDIVE_TSV}" "${WORK_DIR}/bacdive_range.spec" ${Q_BACDIVE_RANGE[@]+"${Q_BACDIVE_RANGE[@]}"}
    write_spec match "${BACDIVE_TSV}" "${WORK_DIR}/bacdive_class.spec" ${Q_BACDIVE_CLASS[@]+"${Q_BACDIVE_CLASS[@]}"}
    for FILTER in ${Q_BACDIVE_RANGE[@]+"${Q_BACDIVE_RANGE[@]}"}; do
        add_filter_name "$(bacdive_filter_name range "${FILTER%%=*}")" "BacDive ${FILTER%%=*} ${FILTER#*=}"
    done
    for FILTER in ${Q_BACDIVE_CLASS[@]+"${Q_BACDIVE_CLASS[@]}"}; do
        add_filter_name "$(bacdive_filter_name class "${FILTER%%=*}")" "BacDive ${FILTER%%=*} ${FILTER#*=}"
    done

    BD_LINKS_FILE="${WORK_DIR}/bacdive_links.tsv"
    awk -F'\t' -v OFS='\t' -v c_taxid="${BD_COL_TAXID}" -v c_rank="${BD_COL_RANK}" \
        -v range_spec="${WORK_DIR}/bacdive_range.spec" -v class_spec="${WORK_DIR}/bacdive_class.spec" \
        -v stats="${WORK_DIR}/stage1.stats" -v summary="${WORK_DIR}/stage1.summary" "${AWK_LIB}"'
    BEGIN {
        nr = load_spec(range_spec, r_col, r_lo, r_hi, r_label)
        nc = load_spec(class_spec, c_col, c_classes, c_unused, c_label)
        for (i = 1; i <= nc; i++) {
            m = split(c_classes[i], want, ",")
            for (k = 1; k <= m; k++) class_want[i, want[k]] = 1
        }
    }
    # cell "coccus=0.25;rod=0.75": the classes recorded for the TaxID
    function class_match(cell, i,   parts, n, j, tok) {
        if (is_empty(cell)) return 0
        n = split(tolower(cell), parts, ";")
        for (j = 1; j <= n; j++) {
            tok = parts[j]
            sub(/=.*/, "", tok)
            gsub(/^[ \t]+|[ \t]+$/, "", tok)
            if ((i, tok) in class_want) return 1
        }
        return 0
    }
    NR == 1 { next }
    {
        if (is_empty($c_taxid)) next
        n_in++
        rank = (is_empty($c_rank) ? "NULL" : $c_rank)
        flags = ""; alive = 1
        for (i = 1; i <= nr; i++) { ok = any_in_range($(r_col[i]), r_lo[i], r_hi[i]); flags = flags ok; alive = alive && ok }
        for (i = 1; i <= nc; i++) { ok = class_match($(c_col[i]), i); flags = flags ok; alive = alive && ok }
        if (above_species(rank)) { n_high++; if (alive) high_pass[rank]++; next }
        n_linkable++
        for (i = 1; i <= nr + nc; i++) { if (substr(flags, i, 1) != "1") break; n_f[i]++ }
        if (alive) pass_rank[rank]++
        print $c_taxid, (rank == "species" ? "species" : "exact"), rank, flags
    }
    END {
        printf "  NCBI TaxIDs in the BacDive table\t%d\n", n_in > stats
        if (n_high) printf "  ... above the species rank (not linked to genomes)\t%d\n", n_high > stats
        printf "  ... species or below (linked to genomes)\t%d\n", n_linkable > stats
        for (i = 1; i <= nr; i++) printf "  ... after %s\t%d\n", r_label[i], n_f[i] > stats
        for (i = 1; i <= nc; i++) printf "  ... after %s\t%d\n", c_label[i], n_f[nr + i] > stats
        # summary: "pass <TAB> <rank> <TAB> <count>" and "high <TAB> <rank> <TAB> <count>"
        for (r in pass_rank) printf "pass\t%s\t%d\n", r, pass_rank[r] > summary
        for (r in high_pass) printf "high\t%s\t%d\n", r, high_pass[r] > summary
        close(summary)
    }' "${BACDIVE_TSV}" > "${BD_LINKS_FILE}"
    touch "${WORK_DIR}/stage1.summary"

    log_stats "${WORK_DIR}/stage1.stats"
    # e.g. "231 NCBI TaxIDs match all the BacDive filters: 200 species, 31 strain"
    N_BD_MATCH="$(awk -F'\t' '$1 == "pass" { n += $3 } END { print n + 0 }' "${WORK_DIR}/stage1.summary")"
    BY_RANK="$(awk -F'\t' '$1 == "pass" { print $3 "\t" $2 }' "${WORK_DIR}/stage1.summary" | sort -k1,1nr | awk -F'\t' '{ printf "%s%d %s", (NR > 1 ? ", " : ""), $1, $2 }')"
    log "message" "Stage 1 complete: ${N_BD_MATCH} NCBI TaxIDs match all the BacDive filters${BY_RANK:+: ${BY_RANK}}."
    HIGH="$(awk -F'\t' '$1 == "high" { print $3 "\t" $2 }' "${WORK_DIR}/stage1.summary" | sort -k1,1nr | awk -F'\t' '{ printf "%s%d %s", (NR > 1 ? ", " : ""), $1, $2 }')"
    if [[ -n "${HIGH}" ]]; then
        log "message" "  TaxIDs above the species rank that match the filters but are not linked to genomes: ${HIGH}."
    fi
    if [[ "${N_BD_MATCH}" -eq 0 ]]; then
        log "warning" "No BacDive TaxID matched all the phenotypic filters, so no genome can be selected."
    fi
else
    log "message" "Stage 1 skipped: no BacDive filter."
fi

############ STAGE 2: GTDB ################
#   GTDB_TAXON: the genome belongs to one of the taxa (any rank). GTDB_RANGE/GTDB_MATCH: GTDB metadata.
#   The keys are accession numbers (GCA_/GCF_ prefix and version removed), because GTDB lists each genome under only
#   one of its GenBank/RefSeq accessions.

GTDB_TAXON_FILE=""   # <key> <TAB> taxa of the query found in its GTDB lineage
GTDB_META_FILE=""    # <key> <TAB> one 0/1 flag per GTDB_RANGE/GTDB_MATCH filter
if [[ "${HAS_GTDB_FILTER}" -eq 1 ]]; then
    log "message" "Stage 2: applying the GTDB filters..."
    GTDB_DIR="${REF_DIR}/GTDB"
    GTDB_TAX_TABLES=()
    GTDB_META_TABLES=()
    for DOMAIN in "${PF_DOMAINS[@]}"; do
        GTDB_TAX_TABLES+=("$(pf_gtdb_taxonomy_table "${GTDB_DIR}" "${DOMAIN}")")
        GTDB_META_TABLES+=("$(pf_gtdb_metadata_table "${GTDB_DIR}" "${DOMAIN}")")
    done

    if [[ "${#Q_GTDB_TAXON[@]}" -gt 0 ]]; then
        for TABLE in "${GTDB_TAX_TABLES[@]}"; do
            [[ -f "${TABLE}" ]] || die "GTDB taxonomy table not found: ${TABLE}. Was the GTDB reference database downloaded?"
        done
        add_filter_name "GTDB_TAXON" "GTDB taxon (${Q_GTDB_TAXON[*]})"
        GTDB_TAXON_FILE="${WORK_DIR}/gtdb_taxon_keys.tsv"
        printf '%s\n' "${Q_GTDB_TAXON[@]}" > "${WORK_DIR}/gtdb_taxa.txt"
        awk -F'\t' -v OFS='\t' -v taxa_file="${WORK_DIR}/gtdb_taxa.txt" "${AWK_LIB}"'
            BEGIN {
                while ((getline t < taxa_file) > 0) { v = tolower(t); sub(/^[a-z]__/, "", v); want[v] = (v in want ? want[v] ";" t : t) }
                close(taxa_file)
            }
            FNR == 1 { delete rank_col; for (i = 1; i <= NF; i++) if ($i ~ /^GTDB_(DOMAIN|PHYLUM|CLASS|ORDER|FAMILY|GENUS|SPECIES)$/) rank_col[i] = 1; next }
            {
                found = ""
                for (i in rank_col) if (tolower($i) in want) found = found (found == "" ? "" : ";") want[tolower($i)]
                if (found != "") print acc_base($1), found
            }
        ' "${GTDB_TAX_TABLES[@]}" > "${GTDB_TAXON_FILE}"
        log "message" "  GTDB genomes in the requested taxa (${Q_GTDB_TAXON[*]}): $(count_lines "${GTDB_TAXON_FILE}")"
    fi

    if [[ $(( ${#Q_GTDB_RANGE[@]} + ${#Q_GTDB_MATCH[@]} )) -gt 0 ]]; then
        for TABLE in "${GTDB_META_TABLES[@]}"; do
            [[ -f "${TABLE}" ]] || die "GTDB metadata table not found: ${TABLE}. Was the GTDB reference database downloaded?"
        done
        write_spec range "${GTDB_META_TABLES[0]}" "${WORK_DIR}/gtdb_range.spec" ${Q_GTDB_RANGE[@]+"${Q_GTDB_RANGE[@]}"}
        write_spec match "${GTDB_META_TABLES[0]}" "${WORK_DIR}/gtdb_match.spec" ${Q_GTDB_MATCH[@]+"${Q_GTDB_MATCH[@]}"}
        for FILTER in ${Q_GTDB_RANGE[@]+"${Q_GTDB_RANGE[@]}"}; do
            add_filter_name "GTDB_${FILTER%%=*}" "GTDB ${FILTER%%=*} ${FILTER#*=}"
        done
        for FILTER in ${Q_GTDB_MATCH[@]+"${Q_GTDB_MATCH[@]}"}; do
            add_filter_name "GTDB_MATCH_${FILTER%%=*}" "GTDB ${FILTER}"
        done
        GTDB_META_FILE="${WORK_DIR}/gtdb_metadata_keys.tsv"
        awk -F'\t' -v OFS='\t' -v range_spec="${WORK_DIR}/gtdb_range.spec" -v match_spec="${WORK_DIR}/gtdb_match.spec" \
            -v stats="${WORK_DIR}/stage2.stats" "${AWK_LIB}"'
            BEGIN {
                nr = load_spec(range_spec, r_col, r_lo, r_hi, r_label)
                nm = load_spec(match_spec, m_col, m_text, m_unused, m_label)
            }
            FNR == 1 { next }
            {
                n_in++
                flags = ""; any = 0
                for (i = 1; i <= nr; i++) { ok = any_in_range($(r_col[i]), r_lo[i], r_hi[i]); flags = flags ok; any = any || ok }
                for (i = 1; i <= nm; i++) { ok = !is_empty($(m_col[i])) && index(tolower($(m_col[i])), m_text[i]) > 0; flags = flags ok; any = any || ok }
                for (i = 1; i <= nr + nm; i++) { if (substr(flags, i, 1) != "1") break; n_f[i]++ }
                if (any) print acc_base($1), flags
            }
            END {
                printf "  GTDB genomes (metadata)\t%d\n", n_in > stats
                for (i = 1; i <= nr; i++) printf "  ... after %s\t%d\n", r_label[i], n_f[i] > stats
                for (i = 1; i <= nm; i++) printf "  ... after %s\t%d\n", m_label[i], n_f[nr + i] > stats
            }
        ' "${GTDB_META_TABLES[@]}" > "${GTDB_META_FILE}"
        log_stats "${WORK_DIR}/stage2.stats"
    fi
    log "message" "Stage 2 complete."
else
    log "message" "Stage 2 skipped: no GTDB filter."
fi

############ STAGE 3: NCBI ################
#   Every assembly is checked against every filter (the BacDive and GTDB filters through stages 1 and 2). The
#   assemblies that satisfy all of them are selected; those that satisfy at least one go to the filters table.

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

# Taxon names (TAXON_NAME) -> TaxIDs: <taxid> <TAB> <name>
TAXON_NAME_LIST=""
if [[ "${#Q_TAXON_NAME[@]}" -gt 0 ]]; then
    require_programs taxonkit
    TAXDUMP_IN_USE="$(pf_ncbi_taxdump_dir "${NCBI_DIR}")"
    [[ -f "${TAXDUMP_IN_USE}/names.dmp" ]] || TAXDUMP_IN_USE=""
    TAXON_NAME_LIST="${WORK_DIR}/taxon_names.tsv"
    resolve_taxon_names "${TAXDUMP_IN_USE}" "${TAXON_NAME_LIST}"
    log "message" "  Taxon names converted to NCBI TaxIDs: $(awk -F'\t' '{ printf "%s%s = %s", (NR > 1 ? ", " : ""), $2, $1 }' "${TAXON_NAME_LIST}")"
fi

# Names of the NCBI filters, in the order in which they are checked below
[[ -n "${ACC_LIST}" ]]          && add_filter_name "ACCESSION" "accession filter"
[[ -n "${TAXID_LIST}" ]]        && add_filter_name "TAXID" "TaxID filter"
[[ -n "${TAXON_LIST}" ]]        && add_filter_name "NCBI_TAXON" "NCBI taxon filter"
[[ -n "${TAXON_NAME_LIST}" ]]   && add_filter_name "NCBI_TAXON_LIST" "NCBI taxon list filter"
[[ -n "${ISO_LIST}" ]]          && add_filter_name "ISOLATION_SOURCE" "isolation source filter"
[[ -n "${COUNTRY_LIST}" ]]      && add_filter_name "COUNTRY" "country filter"
[[ "${Q_GENOME_TYPE}" != "all" ]] && add_filter_name "GENOME_TYPE" "genome type filter"
for FILTER in ${Q_NCBI_RANGE[@]+"${Q_NCBI_RANGE[@]}"}; do
    if [[ "${FILTER%%=*}" == "ASSEMBLY_STATS_GC_PERCENT" ]]; then
        add_filter_name "GC" "NCBI ${FILTER%%=*} ${FILTER#*=}"
    else
        add_filter_name "NCBI_${FILTER%%=*}" "NCBI ${FILTER%%=*} ${FILTER#*=}"
    fi
done
for FILTER in ${Q_NCBI_MATCH[@]+"${Q_NCBI_MATCH[@]}"}; do
    add_filter_name "NCBI_MATCH_${FILTER%%=*}" "NCBI ${FILTER}"
done

SELECTED_TMP="${WORK_DIR}/selected.tsv"
FILTERS_TMP="${WORK_DIR}/filters.tsv"

awk -F'\t' -v OFS='\t' \
    -v c_acc="${C_ACC}" -v c_org="${C_ORG}" -v c_taxid="${C_TAXID}" -v c_iso="${C_ISO}" -v c_country="${C_COUNTRY}" \
    -v c_gtype="${C_GTYPE}" -v c_species_taxid="${C_SPECIES_TAXID}" -v c_strain_taxid="${C_STRAIN_TAXID}" \
    -v name_cols="${LINEAGE_NAME_COLS}" -v taxid_cols="${LINEAGE_TAXID_COLS}" \
    -v bd_file="${BD_LINKS_FILE}" -v gtdb_taxon_file="${GTDB_TAXON_FILE}" -v gtdb_meta_file="${GTDB_META_FILE}" \
    -v acc_file="${ACC_LIST}" -v taxid_file="${TAXID_LIST}" -v taxon_file="${TAXON_LIST}" \
    -v tname_file="${TAXON_NAME_LIST}" -v iso_file="${ISO_LIST}" -v country_file="${COUNTRY_LIST}" \
    -v range_spec="${WORK_DIR}/ncbi_range.spec" -v match_spec="${WORK_DIR}/ncbi_match.spec" \
    -v names_file="${FILTER_NAMES_FILE}" -v filters_out="${FILTERS_TMP}" \
    -v gtype="${Q_GENOME_TYPE}" -v stats="${WORK_DIR}/stage3.stats" "${AWK_LIB}"'
function norm(s) { s = tolower(s); gsub(/[^a-z0-9]+/, " ", s); sub(/ +$/, "", s); return " " s }
# list[i] = the value to search (lower case, or normalized words); orig[i] = the value as given by the user
function load_list(file, list, orig, as_words,   line, n) {
    if (file == "") return 0
    while ((getline line < file) > 0) if (line != "") { n++; list[n] = (as_words ? norm(line) : tolower(line)); orig[n] = line }
    close(file)
    return n + 0
}
# values of the list contained in <text>, joined with ";" (empty if none)
function matches(text, list, orig, n,   i, out) {
    out = ""
    for (i = 1; i <= n; i++) if (index(text, list[i]) > 0) out = out (out == "" ? "" : ";") orig[i]
    return out
}
function add_unique(list, v) { return (index(";" list ";", ";" v ";") > 0 ? list : list (list == "" ? "" : ";") v) }
function or_null(x) { return (x == "" ? "NULL" : x) }
BEGIN {
    has_bd = (bd_file != ""); has_gtdb_taxon = (gtdb_taxon_file != ""); has_gtdb_meta = (gtdb_meta_file != "")
    has_acc = (acc_file != ""); has_taxid = (taxid_file != ""); has_taxon = (taxon_file != "")
    has_tname = (tname_file != ""); has_iso = (iso_file != ""); has_country = (country_file != "")
    has_gtype = (gtype != "all")

    # BacDive TaxIDs: rank species -> linked by the species TaxID; below species -> by the organism/strain TaxID
    if (has_bd) {
        while ((getline line < bd_file) > 0) {
            split(line, f, "\t")
            if (f[2] == "species") bd_species[f[1]] = f[4]; else bd_exact[f[1]] = f[4]
        }
        close(bd_file)
    }
    if (has_gtdb_taxon) { while ((getline line < gtdb_taxon_file) > 0) { split(line, f, "\t"); gtdb_taxon[f[1]] = f[2] } close(gtdb_taxon_file) }
    if (has_gtdb_meta)  { while ((getline line < gtdb_meta_file) > 0)  { split(line, f, "\t"); gtdb_meta[f[1]] = f[2] } close(gtdb_meta_file) }

    load_set(taxid_file, wanted_taxids)
    if (has_acc) { while ((getline a < acc_file) > 0) wanted_acc[acc_number(a)] = 1; close(acc_file) }
    if (has_taxon) { while ((getline a < taxon_file) > 0) wanted_taxon[tolower(a)] = a; close(taxon_file) }
    if (has_tname) { while ((getline line < tname_file) > 0) { split(line, f, "\t"); tname[f[1]] = (f[1] in tname ? tname[f[1]] ";" f[2] : f[2]) } close(tname_file) }
    n_iso = load_list(iso_file, iso_list, iso_orig, 0)
    n_country = load_list(country_file, country_list, country_orig, 1)
    nr = load_spec(range_spec, r_col, r_lo, r_hi, r_label)
    nm = load_spec(match_spec, m_col, m_text, m_unused, m_label)
    n_name_cols = split(name_cols, name_col, ",")
    n_taxid_cols = split(taxid_cols, taxid_col, ",")

    # SELECT_BY names and labels of the filters, in the order checked below
    while ((getline line < names_file) > 0) { split(line, f, "\t"); nf++; fname[nf] = f[1]; flabel[nf] = f[2] }
    close(names_file)
    # Number of filters of stages 1 and 2 (counted from their labels, also when no TaxID/genome passed them)
    for (k = 1; k <= nf; k++) {
        if (flabel[k] ~ /^BacDive /) n_bd++
        else if (flabel[k] ~ /^GTDB / && flabel[k] !~ /^GTDB taxon /) n_gm++
    }

    # Columns with the value that matched each filter with several values
    if (has_taxid)      mcol[++nmc] = "TAXID_MATCH"
    if (has_taxon)      mcol[++nmc] = "NCBI_TAXON_MATCH"
    if (has_tname)      mcol[++nmc] = "NCBI_TAXON_LIST_MATCH"
    if (has_gtdb_taxon) mcol[++nmc] = "GTDB_TAXON_MATCH"
    if (has_iso)        mcol[++nmc] = "ISOLATION_SOURCE_KEYWORD_MATCH"
    if (has_country)    mcol[++nmc] = "COUNTRY_KEYWORD_MATCH"
}
FNR == 1 {
    if (NR == 1) {
        mh = ""
        for (j = 1; j <= nmc; j++) mh = mh OFS mcol[j]
        print $0 mh
        out = "ASSEMBLY_ACCESSION" OFS "ORGANISM_NAME" OFS "SELECTED"
        for (k = 1; k <= nf; k++) out = out OFS fname[k]
        print out mh > filters_out
    }
    next
}
{
    n_in++
    k = 0
    delete mval

    # 1. BacDive (stage 1): the TaxID of the organism/strain (rank below species) or of the species of the assembly
    if (has_bd) {
        flags = ""
        if ($c_taxid in bd_exact) flags = bd_exact[$c_taxid]
        else if ($c_strain_taxid in bd_exact) flags = bd_exact[$c_strain_taxid]
        else if ($c_species_taxid in bd_species) flags = bd_species[$c_species_taxid]
        for (j = 1; j <= n_bd; j++) pass[++k] = (substr(flags, j, 1) == "1")
    }

    # 2. GTDB (stage 2)
    key = acc_base($c_acc)
    if (has_gtdb_taxon) { pass[++k] = (key in gtdb_taxon); mval["GTDB_TAXON_MATCH"] = (key in gtdb_taxon ? gtdb_taxon[key] : "") }
    if (n_gm > 0) {
        flags = (key in gtdb_meta ? gtdb_meta[key] : "")
        for (j = 1; j <= n_gm; j++) pass[++k] = (substr(flags, j, 1) == "1")
    }

    # 3. Accessions (version optional; GenBank/RefSeq pairs match each other)
    if (has_acc) pass[++k] = ((acc_number($c_acc) in wanted_acc) || (acc_base($c_acc) in wanted_acc))

    # 4. TaxID of the organism or of any rank of its lineage
    if (has_taxid) {
        m = ""
        for (i = 1; i <= n_taxid_cols; i++) if ($(taxid_col[i]) in wanted_taxids) m = add_unique(m, $(taxid_col[i]))
        pass[++k] = (m != ""); mval["TAXID_MATCH"] = m
    }

    # 5. NCBI taxon names (lineage names) and 6. list of taxon names (converted to TaxIDs)
    if (has_taxon) {
        m = ""
        for (i = 1; i <= n_name_cols; i++) if (tolower($(name_col[i])) in wanted_taxon) m = add_unique(m, wanted_taxon[tolower($(name_col[i]))])
        pass[++k] = (m != ""); mval["NCBI_TAXON_MATCH"] = m
    }
    if (has_tname) {
        m = ""
        for (i = 1; i <= n_taxid_cols; i++) if ($(taxid_col[i]) in tname) m = add_unique(m, tname[$(taxid_col[i])])
        pass[++k] = (m != ""); mval["NCBI_TAXON_LIST_MATCH"] = m
    }

    # 7. Isolation source (substring) and 8. country (start of a word)
    if (has_iso) {
        m = (is_empty($c_iso) ? "" : matches(tolower($c_iso), iso_list, iso_orig, n_iso))
        pass[++k] = (m != ""); mval["ISOLATION_SOURCE_KEYWORD_MATCH"] = m
    }
    if (has_country) {
        m = (is_empty($c_country) ? "" : matches(norm($c_country), country_list, country_orig, n_country))
        pass[++k] = (m != ""); mval["COUNTRY_KEYWORD_MATCH"] = m
    }

    # 9. Genome type
    if (has_gtype) pass[++k] = (tolower($c_gtype) == gtype)

    # 10. Ranges and text matches on any NCBI column
    for (i = 1; i <= nr; i++) pass[++k] = any_in_range($(r_col[i]), r_lo[i], r_hi[i])
    for (i = 1; i <= nm; i++) pass[++k] = (!is_empty($(m_col[i])) && index(tolower($(m_col[i])), m_text[i]) > 0)

    # Counts after each filter (in order) and the result of the assembly
    any = 0; all = 1
    for (j = 1; j <= k; j++) {
        if (pass[j]) any = 1; else all = 0
        if (all) n_after[j]++
    }
    if (!any) next

    mv = ""
    for (j = 1; j <= nmc; j++) mv = mv OFS or_null(mval[mcol[j]])
    out = $c_acc OFS $c_org OFS (all ? "TRUE" : "FALSE")
    for (j = 1; j <= k; j++) out = out OFS (pass[j] ? "TRUE" : "FALSE")
    print out mv > filters_out
    if (all) print $0 mv
}
END {
    printf "  NCBI genome assemblies searched\t%d\n", n_in > stats
    for (j = 1; j <= nf; j++) printf "  ... after %s\t%d\n", flabel[j], n_after[j] > stats
}' "${NCBI_TABLES[@]}" > "${SELECTED_TMP}"

log_stats "${WORK_DIR}/stage3.stats"
TOTAL_SELECTED="$(count_records "${SELECTED_TMP}")"
log "message" "Stage 3 complete: ${TOTAL_SELECTED} genome assemblies selected; $(count_records "${FILTERS_TMP}") satisfy at least one filter."

############ BACDIVE ANNOTATION ################
#   The BacDive traits of each selected assembly are added as BACDIVE_* columns, from the TaxID table (strains of the
#   same TaxID merged). An assembly uses the row of its organism or strain TaxID when it is below the species rank,
#   otherwise the row of its species TaxID. BACDIVE_TAXID and BACDIVE_LINK_LEVEL (rank) tell which row was used.

if [[ -f "${BACDIVE_TSV}" && "${TOTAL_SELECTED}" -gt 0 ]]; then
    log "message" "Adding the BacDive traits of the selected genomes..."
    BD_SPEC=""
    for TRAIT in bacdive_strains ${BACDIVE_ANNOTATION_TRAITS[@]+"${BACDIVE_ANNOTATION_TRAITS[@]}"}; do
        INDEX="$(get_column_index "${BACDIVE_TSV}" "${TRAIT}")"
        if [[ -n "${INDEX}" ]]; then
            BD_SPEC="${BD_SPEC:+${BD_SPEC},}${INDEX}=${TRAIT}"
        fi
    done

    awk -F'\t' -v OFS='\t' -v bd_file="${BACDIVE_TSV}" \
        -v bd_taxid_col="$(require_column "${BACDIVE_TSV}" ncbi_tax_id)" -v bd_rank_col="$(require_column "${BACDIVE_TSV}" ncbi_tax_id_rank)" \
        -v spec="${BD_SPEC}" -v c_acc="${C_ACC}" -v c_taxid="${C_TAXID}" -v c_strain_taxid="${C_STRAIN_TAXID}" \
        -v c_species_taxid="${C_SPECIES_TAXID}" "${AWK_LIB}"'
    BEGIN {
        ns = split(spec, sp, ",")
        for (k = 1; k <= ns; k++) { split(sp[k], p, "="); col[k] = p[1] + 0; name[k] = (p[2] == "bacdive_strains" ? "BACDIVE_STRAINS" : "BACDIVE_" toupper(p[2])) }
        while ((getline line < bd_file) > 0) {
            if (++rows == 1) continue
            split(line, f, "\t")
            t = f[bd_taxid_col]; r = f[bd_rank_col]
            if (is_empty(t) || above_species(r)) continue
            v = ""
            for (k = 1; k <= ns; k++) v = v OFS (f[col[k]] == "" ? "NULL" : f[col[k]])
            if (r == "species") species_row[t] = v; else exact_row[t] = v
            rank[t] = r
        }
        close(bd_file)
        out = "ASSEMBLY_ACCESSION" OFS "BACDIVE_TAXID" OFS "BACDIVE_LINK_LEVEL"
        for (k = 1; k <= ns; k++) out = out OFS name[k]
        print out
        empty = ""
        for (k = 1; k <= ns; k++) empty = empty OFS "NULL"
    }
    NR == 1 { next }
    {
        if ($c_taxid in exact_row)                t = $c_taxid
        else if ($c_strain_taxid in exact_row)    t = $c_strain_taxid
        else if ($c_species_taxid in species_row) t = $c_species_taxid
        else                                      t = ""
        if (t == "") print $c_acc, "NULL", "NULL" empty
        else print $c_acc, t, rank[t] ((t in exact_row) ? exact_row[t] : species_row[t])
    }' "${SELECTED_TMP}" > "${WORK_DIR}/bacdive_annotation.tsv"

    add_metadata_columns "${SELECTED_TMP}" "${WORK_DIR}/bacdive_annotation.tsv"
    log "message" "  Selected genomes linked to BacDive: $(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "BACDIVE_LINK_LEVEL") c = i; next } $c != "NULL" { n[$c]++; t++ } END { printf "%d", t; s = ""; for (r in n) s = s (s == "" ? "" : ", ") n[r] " " r; if (s != "") printf " (%s)", s }' "${SELECTED_TMP}")"
fi

mkdir -p "$(dirname "${OUTPUT_TSV}")" "$(dirname "${FILTERS_TSV}")"
mv "${SELECTED_TMP}" "${OUTPUT_TSV}"
mv "${FILTERS_TMP}" "${FILTERS_TSV}"
log "message" "Selected genome metadata (one row per assembly) saved to: ${OUTPUT_TSV}"
log "message" "Assemblies that satisfy at least one filter (TRUE/FALSE per filter) saved to: ${FILTERS_TSV}"
