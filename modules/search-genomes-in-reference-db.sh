#!/usr/bin/env bash
#
# search-genomes-in-reference-db.sh - Prokfinder module 1: genome search in the reference database
#
# Reads a query file (genome-selection-query.txt, written by "main-prokfinder.sh search-genomes") and selects, in
# the NCBI tables, the genome assemblies that satisfy the selection expression. Each filter is one option of the
# command line (several values of the same option = OR); the filters are combined by the expression, with AND, OR
# and parentheses (AND is evaluated before OR). The filters can use the metadata of NCBI, the phenotypic traits of
# BacDive and the taxonomy of GTDB:
#
#   BacDive -> the BacDive table used has one row per NCBI TaxID (the strains of the same TaxID merged, see
#              download-reference-db-bacdive.sh). A TaxID of rank species is linked to the assemblies of that species
#              (NCBI_SPECIES_TAXID); a TaxID below the species rank (strain, subspecies...) to the assemblies of that
#              organism or strain (exact TaxID); TaxIDs above the species rank (genus, family...) are not linked. An
#              assembly linked to a TaxID below the species rank uses its row; otherwise the row of its species.
#   GTDB    -> the GTDB lineage of the assembly (GTDB lists each genome under only one of its GenBank/RefSeq
#              accessions, so the accession number without prefix is used)
#   NCBI    -> every assembly of the NCBI table (GenBank by default) is checked against every filter, and the
#              expression decides whether it is selected
#
# QUERY FILE (TAB-separated; empty lines and lines starting with "#" are ignored):
#   ASSEMBLY_SOURCE <TAB> genbank|refseq                 NCBI table that is searched (default: genbank)
#   FILTER <TAB> <id> <TAB> <KEY> <TAB> <value>          one line per value; lines with the same id are one filter
#   EXPRESSION <TAB> <id> AND ( <id> OR <id> ) ...       how the filters are combined (default: AND of all filters)
#
#   KEY (values of the same filter combined with OR)
#     ACCESSION         GCA_/GCF_ accession with version (GCA_x.1 also matches its paired GCF_x.1)
#     TAXID             NCBI TaxID of the organism or of any rank of its lineage (descendant taxa included)
#     NCBI_TAXON        name of a taxon of the NCBI lineage (any rank, exact match, case-insensitive)
#     TAXON_NAME        name of an NCBI taxon, converted to its TaxID with "taxonkit name2taxid" (Bacteria and
#                       Archaea only); then the same as TAXID. A name without TaxID is an error
#     ISOLATION_SOURCE  extended regular expression (POSIX ERE; a plain text is one too) found in the isolation
#                       source (case-insensitive)
#     COUNTRY           country/place at the start of a word of the geographic location (case-insensitive)
#     GTDB_TAXON        name of a taxon of the GTDB lineage (any rank; "g__Thermus" or "Thermus")
#     BACDIVE_CLASS     <trait>=<class>[,<class>...]: at least one of the classes was recorded for the TaxID
#   KEY (single value)
#     GENOME_TYPE       isolate_genome or mag
#     GC                <min>:<max> GC content (%) of the assembly (either side may be empty)
#     BACDIVE_RANGE     <column>=<min>:<max>: the value of the TaxID (mean of its strains) is inside the range
#
# Usage:
#   search-genomes-in-reference-db.sh -q <query_file> --validate-only [--taxdump-dir <dir>]
#   search-genomes-in-reference-db.sh -q <query_file> -r <reference_db_dir> -o <selected_metadata.tsv>
#
# Output: <selected_metadata.tsv>, one row per selected assembly: NCBI columns, one SELECTED_BY_* column (TRUE/FALSE)
# per filter, the *_MATCH columns (values that matched each filter with several values, e.g.
# ISOLATION_SOURCE_KEYWORD_MATCH) and the BACDIVE_* columns.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

usage() {
    cat << EOF
Usage: $(basename "$0") -q <query_file> [-r <reference_db_dir> -o <output_tsv>] [--validate-only]

  -q <file>             Query file (see the header of this script for the format)
  -r <dir>              Reference databases folder (~/prokfinder-refdb)
  -o <file>             Output table with the selected genome assemblies
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
TAXDUMP_DIR=""
VALIDATE_ONLY=0

while [[ $# -gt 0 ]]; do
    case "$1" in
        -q) QUERY_FILE="$2"; shift 2 ;;
        -r) REF_DIR="${2%/}"; shift 2 ;;
        -o) OUTPUT_TSV="$2"; shift 2 ;;
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

Q_ASSEMBLY_SOURCE="genbank"
EXPRESSION=""
F_IDS=()                 # filter ids, in the order of the query file
declare -A F_KEY=()      # id -> key
declare -A F_VALUES=()   # id -> values, one per line
declare -A F_SEEN=()     # id TAB value -> 1
BAD_ACCESSIONS=()

# Keys whose filter has a single value
SINGLE_VALUE_KEYS=" GENOME_TYPE GC BACDIVE_RANGE BACDIVE_CLASS "

LINE_NUMBER=0
while IFS= read -r LINE || [[ -n "${LINE}" ]]; do
    LINE_NUMBER=$((LINE_NUMBER + 1))
    LINE="${LINE//$'\r'/}"
    if [[ -z "${LINE//[[:space:]]/}" || "${LINE}" =~ ^[[:space:]]*# ]]; then
        continue
    fi
    IFS=$'\t' read -r KIND ID KEY VALUE <<< "${LINE}"
    case "${KIND}" in
        ASSEMBLY_SOURCE)
            Q_ASSEMBLY_SOURCE="${ID,,}"
            if [[ "${Q_ASSEMBLY_SOURCE}" != "genbank" && "${Q_ASSEMBLY_SOURCE}" != "refseq" ]]; then
                die "Invalid assembly source '${ID}'. Must be 'genbank' or 'refseq'."
            fi
            continue
            ;;
        EXPRESSION)
            EXPRESSION="${LINE#*$'\t'}"
            continue
            ;;
        FILTER) ;;
        *) die "Unknown line ${LINE_NUMBER} of the query file: '${LINE}'" ;;
    esac

    if [[ ! "${ID}" =~ ^[A-Za-z0-9_]+$ ]]; then
        die "Query line ${LINE_NUMBER}: invalid filter id '${ID}'."
    fi
    VALUE="${VALUE:-}"
    VALUE="${VALUE#"${VALUE%%[![:space:]]*}"}"
    VALUE="${VALUE%"${VALUE##*[![:space:]]}"}"
    if [[ -z "${VALUE}" ]]; then
        die "Query line ${LINE_NUMBER} (${KEY}) has no value."
    fi
    if [[ -n "${F_KEY[${ID}]:-}" && "${F_KEY[${ID}]}" != "${KEY}" ]]; then
        die "Query line ${LINE_NUMBER}: filter ${ID} has two different keys (${F_KEY[${ID}]} and ${KEY})."
    fi
    if [[ -n "${F_KEY[${ID}]:-}" && "${SINGLE_VALUE_KEYS}" == *" ${KEY} "* ]]; then
        die "Query line ${LINE_NUMBER}: filter ${ID} (${KEY}) takes a single value."
    fi

    case "${KEY}" in
        ACCESSION)
            VALUE="${VALUE^^}"
            if [[ ! "${VALUE}" =~ ^GC[AF]_[0-9]+[.][0-9]+$ ]]; then
                BAD_ACCESSIONS+=("${VALUE}")
                continue
            fi
            ;;
        TAXID)
            [[ "${VALUE}" =~ ^[0-9]+$ ]] || die "Invalid NCBI TaxID '${VALUE}' (expected a number)."
            ;;
        NCBI_TAXON|TAXON_NAME|ISOLATION_SOURCE|COUNTRY|GTDB_TAXON) ;;
        GENOME_TYPE)
            VALUE="${VALUE,,}"
            [[ "${VALUE}" == "isolate_genome" || "${VALUE}" == "mag" ]] \
                || die "Invalid genome type '${VALUE}'. Must be 'isolate_genome' or 'mag'."
            ;;
        GC)
            validate_range "GC" "${VALUE}"
            ;;
        BACDIVE_RANGE)
            if [[ "${VALUE}" != *=* || -z "${VALUE%%=*}" ]]; then
                die "Invalid BACDIVE_RANGE '${VALUE}'. Expected <column>=<min>:<max>."
            fi
            validate_range "${VALUE%%=*}" "${VALUE#*=}"
            ;;
        BACDIVE_CLASS)
            if [[ "${VALUE}" != *=* ]]; then
                die "Invalid BACDIVE_CLASS '${VALUE}'. Expected <trait>=<class>[,<class>...] (e.g., gram_stain=negative)."
            fi
            TRAIT="$(printf '%s' "${VALUE%%=*}" | tr '[:upper:]' '[:lower:]' | tr -d '[:space:]')"
            if [[ ! "${TRAIT}" =~ ^[a-z0-9_]+$ || -z "${TRAIT_CLASSES[${TRAIT}]:-}" ]]; then
                die "Unknown trait '${VALUE%%=*}'. Traits divided in classes: ${TRAIT_ORDER[*]}."
            fi
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
            VALUE="${TRAIT}=${CLASS_LIST}"
            ;;
        *)
            die "Unknown key '${KEY}' in line ${LINE_NUMBER} of the query file."
            ;;
    esac

    if [[ -z "${F_KEY[${ID}]:-}" ]]; then
        F_IDS+=("${ID}")
        F_KEY["${ID}"]="${KEY}"
        F_VALUES["${ID}"]=""
    fi
    # A value repeated in the same filter (e.g., one line per strain of the same species) is kept once
    [[ -z "${F_SEEN[${ID}$'\t'${VALUE}]:-}" ]] || continue
    F_SEEN["${ID}"$'\t'"${VALUE}"]=1
    F_VALUES["${ID}"]+="${VALUE}"$'\n'
done < "${QUERY_FILE}"

# Every invalid accession is reported at once
if [[ "${#BAD_ACCESSIONS[@]}" -gt 0 ]]; then
    # (awk reads the whole list: "head" would close the pipe early and kill printf with SIGPIPE, exit code 141)
    LIST="$(printf '%s\n' "${BAD_ACCESSIONS[@]}" | awk 'NR <= 20' | paste -sd ' ')"
    MORE=$(( ${#BAD_ACCESSIONS[@]} > 20 ? ${#BAD_ACCESSIONS[@]} - 20 : 0 ))
    die "${#BAD_ACCESSIONS[@]} invalid accession(s): ${LIST}$( (( MORE > 0 )) && echo " (and ${MORE} more)"). An accession is GCA_ or GCF_, the number and the version (e.g., GCA_000005845.2)."
fi
for ID in ${F_IDS[@]+"${F_IDS[@]}"}; do
    [[ -n "${F_VALUES[${ID}]}" ]] || die "Filter ${ID} (${F_KEY[${ID}]}) has no valid value."
done

# The isolation sources are extended regular expressions (a plain text is one too): all of them must compile
ISO_PATTERNS="${WORK_DIR}/isolation_source_patterns.txt"
for ID in ${F_IDS[@]+"${F_IDS[@]}"}; do
    [[ "${F_KEY[${ID}]}" != "ISOLATION_SOURCE" ]] || printf '%s' "${F_VALUES[${ID}]}"
done > "${ISO_PATTERNS}"
if [[ -s "${ISO_PATTERNS}" ]] && ! awk '{ if ("" ~ tolower($0)) n++ }' "${ISO_PATTERNS}" 2> "${WORK_DIR}/regex.err"; then
    die "Invalid regular expression in the isolation sources (-i/-I): $(tr '\n' ' ' < "${WORK_DIR}/regex.err")Characters such as ( ) [ ] . * + ? | ^ \$ have a special meaning; write \\ before them to search them as text (e.g., \"soil \\(forest\\)\")."
fi
if [[ "${#F_IDS[@]}" -eq 0 ]]; then
    die "The query has no selection filter."
fi

############ SELECTION EXPRESSION ################
#   The expression is converted to postfix notation (shunting-yard): AND is evaluated before OR, and parentheses
#   group filters. Every filter must appear in the expression.

if [[ -z "${EXPRESSION//[[:space:]]/}" ]]; then
    EXPRESSION="$(printf '%s\n' "${F_IDS[@]}" | paste -sd ' ' | sed 's/ / AND /g')"
fi
# Usage: expression_error <message>
expression_error() {
    die "Invalid selection expression '${EXPRESSION}': $1"
}
declare -A F_INDEX=()
for (( I = 0; I < ${#F_IDS[@]}; I++ )); do
    F_INDEX["${F_IDS[$I]}"]=$((I + 1))
done
declare -A F_USED=()
RPN=()
STACK=()
PREV="start"   # start | operand | operator | open | close
read -ra TOKENS <<< "${EXPRESSION//(/ ( }"
TOKENS_STR="${TOKENS[*]}"
read -ra TOKENS <<< "${TOKENS_STR//)/ ) }"
for TOKEN in ${TOKENS[@]+"${TOKENS[@]}"}; do
    case "${TOKEN^^}" in
        AND|OR)
            [[ "${PREV}" == "operand" || "${PREV}" == "close" ]] || expression_error "${TOKEN^^} must be between two filters."
            while [[ "${#STACK[@]}" -gt 0 ]]; do
                TOP="${STACK[-1]}"
                if [[ "${TOP}" == "AND" || ( "${TOP}" == "OR" && "${TOKEN^^}" == "OR" ) ]]; then
                    RPN+=("${TOP}"); unset 'STACK[-1]'
                else
                    break
                fi
            done
            STACK+=("${TOKEN^^}")
            PREV="operator"
            ;;
        "(")
            [[ "${PREV}" == "start" || "${PREV}" == "operator" || "${PREV}" == "open" ]] || expression_error "AND or OR is missing before '('."
            STACK+=("(")
            PREV="open"
            ;;
        ")")
            [[ "${PREV}" == "operand" || "${PREV}" == "close" ]] || expression_error "')' must close a group of filters."
            FOUND=0
            while [[ "${#STACK[@]}" -gt 0 ]]; do
                TOP="${STACK[-1]}"; unset 'STACK[-1]'
                if [[ "${TOP}" == "(" ]]; then FOUND=1; break; fi
                RPN+=("${TOP}")
            done
            [[ "${FOUND}" -eq 1 ]] || expression_error "')' without '('."
            PREV="close"
            ;;
        *)
            [[ -n "${F_INDEX[${TOKEN}]:-}" ]] || expression_error "unknown filter '${TOKEN}'."
            [[ "${PREV}" == "start" || "${PREV}" == "operator" || "${PREV}" == "open" ]] || expression_error "AND or OR is missing before ${TOKEN}."
            RPN+=("${F_INDEX[${TOKEN}]}")
            F_USED["${TOKEN}"]=1
            PREV="operand"
            ;;
    esac
done
[[ "${PREV}" == "operand" || "${PREV}" == "close" ]] || expression_error "it must end with a filter or ')'."
while [[ "${#STACK[@]}" -gt 0 ]]; do
    TOP="${STACK[-1]}"; unset 'STACK[-1]'
    [[ "${TOP}" != "(" ]] || expression_error "'(' without ')'."
    RPN+=("${TOP}")
done
for ID in "${F_IDS[@]}"; do
    [[ -n "${F_USED[${ID}]:-}" ]] || expression_error "filter ${ID} is not used."
done

############ TAXON NAMES -> NCBI TAXIDS ################

# Converts every TAXON_NAME name into NCBI TaxIDs of Bacteria or Archaea (taxonkit name2taxid). Stops with the list of
# the names that have no such TaxID. A name with several TaxIDs (homonyms) keeps all of them.
# Usage: resolve_taxon_names <taxdump_dir, or "" for the taxonkit default> <output: taxid TAB name>
resolve_taxon_names() {
    local TAXDUMP="$1"
    local OUT="$2"
    local -a ARGS=()
    local MISSING ID FIRST
    if [[ -n "${TAXDUMP}" ]]; then
        ARGS=(--data-dir "${TAXDUMP}")
    fi

    for ID in "${F_IDS[@]}"; do
        if [[ "${F_KEY[${ID}]}" == "TAXON_NAME" ]]; then
            printf '%s' "${F_VALUES[${ID}]}"
        fi
    done | sort -u > "${WORK_DIR}/taxon_names.txt"

    # name2taxid prints "<name> <TAB> <taxid>" (one line per TaxID; empty TaxID when the name is unknown)
    taxonkit name2taxid ${ARGS[@]+"${ARGS[@]}"} < "${WORK_DIR}/taxon_names.txt" > "${WORK_DIR}/name2taxid.tsv" 2> "${WORK_DIR}/name2taxid.err" \
        || die "taxonkit name2taxid failed: $(head -n 3 "${WORK_DIR}/name2taxid.err" | tr '\n' ' ')"
    awk -F'\t' '$2 ~ /^[0-9]+$/ { print $2 }' "${WORK_DIR}/name2taxid.tsv" | sort -u \
        | taxonkit lineage -t ${ARGS[@]+"${ARGS[@]}"} > "${WORK_DIR}/name2taxid.lineage" 2>> "${WORK_DIR}/name2taxid.err" \
        || die "taxonkit lineage failed: $(head -n 3 "${WORK_DIR}/name2taxid.err" | tr '\n' ' ')"

    # Only TaxIDs whose lineage contains Bacteria (2) or Archaea (2157)
    awk -F'\t' -v OFS='\t' -v lineage="${WORK_DIR}/name2taxid.lineage" '
        BEGIN { while ((getline line < lineage) > 0) { split(line, f, "\t"); if ((";" f[3] ";") ~ /;(2|2157);/) prok[f[1]] = 1 } }
        $2 ~ /^[0-9]+$/ && ($2 in prok) { print $2, $1 }
    ' "${WORK_DIR}/name2taxid.tsv" | sort -u > "${OUT}"

    # The first line of a list file (-N) that is not a taxon is its header (e.g., "Species Name"): it is skipped
    for ID in "${F_IDS[@]}"; do
        [[ "${F_KEY[${ID}]}" == "TAXON_NAME" ]] || continue
        FIRST="${F_VALUES[${ID}]%%$'\n'*}"
        if ! awk -F'\t' -v name="${FIRST,,}" 'tolower($2) == name { found = 1 } END { exit !found }' "${OUT}"; then
            if [[ "${VALIDATE_ONLY}" -eq 0 ]]; then
                log "warning" "First line of the taxon list of filter ${ID} ('${FIRST}') is not a Bacteria/Archaea taxon: taken as a header and skipped."
            fi
            F_VALUES["${ID}"]="${F_VALUES[${ID}]#*$'\n'}"
            [[ -n "${F_VALUES[${ID}]}" ]] || die "Filter ${ID} (taxon list) has no valid value."
        fi
    done

    for ID in "${F_IDS[@]}"; do
        if [[ "${F_KEY[${ID}]}" == "TAXON_NAME" ]]; then
            printf '%s' "${F_VALUES[${ID}]}"
        fi
    done | sort -u > "${WORK_DIR}/taxon_names.txt"
    MISSING="$(awk -F'\t' -v found="${OUT}" '
        BEGIN { while ((getline line < found) > 0) { split(line, f, "\t"); ok[tolower(f[2])] = 1 } }
        !(tolower($0) in ok) && ++n <= 20 { printf "%s\"%s\"", (n > 1 ? ", " : ""), $0 }
        END { if (n > 20) printf " (and %d more)", n - 20 }' "${WORK_DIR}/taxon_names.txt")"
    if [[ -n "${MISSING}" ]]; then
        die "Taxon name(s) not found in the NCBI taxonomy of Bacteria and Archaea: ${MISSING}. Check the spelling (scientific names, e.g., Thermus or Thermus aquaticus)."
    fi
}

HAS_TAXON_NAME=0
for ID in "${F_IDS[@]}"; do
    [[ "${F_KEY[${ID}]}" == "TAXON_NAME" ]] && HAS_TAXON_NAME=1
done

if [[ "${VALIDATE_ONLY}" -eq 1 ]]; then
    # The taxon names are checked now only if taxonkit and a taxdump are available here; otherwise, during the search
    if [[ "${HAS_TAXON_NAME}" -eq 1 ]] && command -v taxonkit >/dev/null 2>&1; then
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

# Taxon names (-N) -> NCBI TaxIDs, now: a header removed from a list must not reach the filter tables built below
if [[ "${HAS_TAXON_NAME}" -eq 1 ]]; then
    require_programs taxonkit
    TAXDUMP_IN_USE="$(pf_ncbi_taxdump_dir "${NCBI_DIR}")"
    [[ -f "${TAXDUMP_IN_USE}/names.dmp" ]] || TAXDUMP_IN_USE=""
    resolve_taxon_names "${TAXDUMP_IN_USE}" "${WORK_DIR}/taxon_names.tsv"
    log "message" "Taxon names converted to NCBI TaxIDs: $(wc -l < "${WORK_DIR}/taxon_names.txt") names, $(cut -f1 "${WORK_DIR}/taxon_names.tsv" | sort -u | wc -l) TaxIDs ($(awk -F'\t' 'NR <= 5 { printf "%s%s = %s", (NR > 1 ? ", " : ""), $2, $1 } END { if (NR > 5) printf ", ..." }' "${WORK_DIR}/taxon_names.tsv"))"
fi

NCBI_TABLES=()
for DOMAIN in "${PF_DOMAINS[@]}"; do
    TABLE="$(pf_ncbi_metadata_table "${NCBI_DIR}" "${DOMAIN}" "${Q_ASSEMBLY_SOURCE}")"
    if [[ ! -f "${TABLE}" ]]; then
        die "NCBI table not found: ${TABLE}. Install the reference databases with: main-prokfinder.sh setup"
    fi
    if [[ "${#NCBI_TABLES[@]}" -gt 0 && "$(head -n 1 "${TABLE}")" != "$(head -n 1 "${NCBI_TABLES[0]}")" ]]; then
        die "The columns of ${TABLE} differ from those of ${NCBI_TABLES[0]}."
    fi
    NCBI_TABLES+=("${TABLE}")
done
NCBI_HEADER_FILE="${NCBI_TABLES[0]}"

# Common awk functions
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
# GCA_000005845.2 -> 000005845.2 (GenBank and RefSeq assemblies of the same genome share the number)
function acc_number(acc,   a) { a = acc; sub(/^GC[AF]_/, "", a); return a }
function acc_base(acc,   a) { a = acc_number(acc); sub(/[.][0-9]+$/, "", a); return a }
# Ranks above species: BacDive TaxIDs of these ranks are not linked to genomes
function above_species(r) {
    return r ~ /^(superkingdom|domain|realm|kingdom|subkingdom|phylum|subphylum|superclass|class|subclass|infraclass|superorder|order|suborder|infraorder|parvorder|superfamily|family|subfamily|tribe|subtribe|genus|subgenus|section|subsection|series|species group|species subgroup)$/
}
'

# Name of a BacDive trait in the SELECTED_BY columns (the option of main-prokfinder.sh)
bacdive_trait_name() {
    case "${1,,}" in
        temperature_opt)        echo "OPT_TEMP" ;;
        ph_opt)                 echo "OPT_PH" ;;
        halophily_opt)          echo "OPT_SALT" ;;
        temperature_preference) echo "TEMP_GROWTH" ;;
        ph_preference)          echo "PH_GROWTH" ;;
        halophily_preference)   echo "SALT_GROWTH" ;;
        cell_shape)             echo "CELL_SHAPE" ;;
        oxygen_tolerance)       echo "OXYGEN" ;;
        *)                      echo "${1^^}" ;;
    esac
}

# Columns written for each filter, in the order of the filters:
#   flags with a single option  SELECTED_BY_<FLAG>-<OPTION>    e.g. SELECTED_BY_GENOME_TYPE-ISOLATE-GENOME
#   numeric flags               SELECTED_BY_<FLAG>-<RANGE>     e.g. SELECTED_BY_OPT_TEMP-45:80
#   flags with several values   SELECTED_BY_<FLAG> + a column with the values that matched (separated by ";")
# A name used twice (the same option given twice) gets a number (_2, _3...).
FILTERS_SPEC="${WORK_DIR}/filters.tsv"     # <index> <key> <arg1> <arg2> <SELECTED_BY name> <match name> <label>
FILTER_VALUES="${WORK_DIR}/filter_values.tsv"   # <index> <value> (filters with a list of values)
: > "${FILTERS_SPEC}"
: > "${FILTER_VALUES}"
declare -A USED_NAME=()
BD_FILTERS=()      # indexes of the BacDive filters, in order
GTDB_FILTERS=()    # indexes of the GTDB filters
for (( I = 0; I < ${#F_IDS[@]}; I++ )); do
    ID="${F_IDS[$I]}"
    IDX=$((I + 1))
    KEY="${F_KEY[${ID}]}"
    VALUES="${F_VALUES[${ID}]%$'\n'}"
    ARG1=""; ARG2=""; MATCH=""
    case "${KEY}" in
        ACCESSION)        NAME="ACCESSION";        LABEL="accession" ;;
        TAXID)            NAME="TAXID";            MATCH="TAXID_MATCH";            LABEL="NCBI TaxID" ;;
        NCBI_TAXON)       NAME="NCBI_TAXON";       MATCH="NCBI_TAXON_MATCH";       LABEL="NCBI taxon" ;;
        TAXON_NAME)       NAME="NCBI_TAXON_LIST";  MATCH="NCBI_TAXON_LIST_MATCH";  LABEL="NCBI taxon list" ;;
        ISOLATION_SOURCE) NAME="ISOLATION_SOURCE"; MATCH="ISOLATION_SOURCE_KEYWORD_MATCH"; LABEL="isolation source" ;;
        COUNTRY)          NAME="COUNTRY";          MATCH="COUNTRY_KEYWORD_MATCH";  LABEL="country" ;;
        GTDB_TAXON)       NAME="GTDB_TAXON";       MATCH="GTDB_TAXON_MATCH";       LABEL="GTDB taxon"; GTDB_FILTERS+=("${IDX}") ;;
        GENOME_TYPE)
            ARG1="${VALUES}"
            NAME="GENOME_TYPE-$(printf '%s' "${VALUES^^}" | tr '_' '-')"; LABEL="genome type" ;;
        GC)
            ARG1="${VALUES%%:*}"; ARG2="${VALUES#*:}"
            NAME="GC-${VALUES}"; LABEL="GC (%)" ;;
        BACDIVE_RANGE)
            ARG1="${VALUES%%=*}"; ARG2="${VALUES#*=}"
            NAME="$(bacdive_trait_name "${ARG1}")-${ARG2}"; LABEL="BacDive ${ARG1}"
            BD_FILTERS+=("${IDX}") ;;
        BACDIVE_CLASS)
            ARG1="${VALUES%%=*}"; ARG2="${VALUES#*=}"
            NAME="$(bacdive_trait_name "${ARG1}")"; MATCH="${NAME}_MATCH"; LABEL="BacDive ${ARG1}"
            BD_FILTERS+=("${IDX}") ;;
    esac
    NAME="SELECTED_BY_${NAME}"
    if [[ -n "${USED_NAME[${NAME}]:-}" ]]; then
        N=2
        while [[ -n "${USED_NAME[${NAME}_${N}]:-}" ]]; do N=$((N + 1)); done
        NAME="${NAME}_${N}"
        [[ -n "${MATCH}" ]] && MATCH="${MATCH}_${N}"
    fi
    USED_NAME["${NAME}"]=1
    case "${KEY}" in
        GENOME_TYPE|GC|BACDIVE_RANGE|BACDIVE_CLASS) LABEL="${LABEL} ${VALUES#*=}" ;;
        *)
            SHOWN="$(printf '%s\n' "${VALUES}" | awk 'NR <= 5 { printf "%s%s", (NR > 1 ? "," : ""), $0 } END { if (NR > 5) printf ",... (%d values)", NR }')"
            LABEL="${LABEL} (${SHOWN})"
            ;;
    esac
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${IDX}" "${KEY}" "${ARG1}" "${ARG2}" "${NAME}" "${MATCH}" "${LABEL}" >> "${FILTERS_SPEC}"
    case "${KEY}" in
        ACCESSION|TAXID|NCBI_TAXON|TAXON_NAME|ISOLATION_SOURCE|COUNTRY|GTDB_TAXON)
            printf '%s\n' "${VALUES}" | awk -v i="${IDX}" -v OFS='\t' 'NF { print i, $0 }' >> "${FILTER_VALUES}" ;;
    esac
done

log "message" "Searching the reference database (NCBI ${Q_ASSEMBLY_SOURCE} tables)..."
log "message" "Selection expression: ${EXPRESSION}"

############ BACDIVE FILTERS ################
#   One row per NCBI TaxID (the BacDive strains of the same TaxID merged). Range filters: the value of the TaxID
#   (mean of its strains) is inside the range. Class filters: at least one of the requested classes was recorded for
#   the TaxID. TaxIDs with no value for a filtered trait do not pass. Output: one line per TaxID that can be linked to
#   genomes: <taxid> <species|exact> <rank> and, per BacDive filter, "1" (range) or the classes found ("" = no).

BD_LINKS_FILE=""
if [[ "${#BD_FILTERS[@]}" -gt 0 ]]; then
    if [[ ! -f "${BACDIVE_TSV}" ]]; then
        die "BacDive filters require the BacDive TaxID table, which was not found: ${BACDIVE_TSV}. Update the reference databases with: main-prokfinder.sh setup"
    fi
    BD_SPEC="${WORK_DIR}/bacdive.spec"   # <filter index> <kind> <column index> <arg1> <arg2> <label>
    : > "${BD_SPEC}"
    for IDX in "${BD_FILTERS[@]}"; do
        # (fields separated by \x1f: unlike TAB, it keeps the empty fields when read)
        IFS=$'\x1f' read -r KEY ARG1 ARG2 LABEL < <(awk -F'\t' -v OFS='\x1f' -v i="${IDX}" '$1 == i { print $2, $3, $4, $7 }' "${FILTERS_SPEC}")
        COL="$(get_column_index "${BACDIVE_TSV}" "${ARG1}")"
        [[ -n "${COL}" ]] || die "Column '${ARG1}' not found in ${BACDIVE_TSV}."
        if [[ "${KEY}" == "BACDIVE_RANGE" ]]; then
            printf '%s\trange\t%s\t%s\t%s\t%s\n' "${IDX}" "${COL}" "${ARG2%%:*}" "${ARG2#*:}" "${LABEL}" >> "${BD_SPEC}"
        else
            printf '%s\tclass\t%s\t%s\t\t%s\n' "${IDX}" "${COL}" "${ARG2}" "${LABEL}" >> "${BD_SPEC}"
        fi
    done

    log "message" "BacDive filters (one row per NCBI TaxID, strains merged):"
    BD_LINKS_FILE="${WORK_DIR}/bacdive_links.tsv"
    awk -F'\t' -v OFS='\t' -v c_taxid="$(require_column "${BACDIVE_TSV}" ncbi_tax_id)" \
        -v c_rank="$(require_column "${BACDIVE_TSV}" ncbi_tax_id_rank)" -v spec="${BD_SPEC}" \
        -v stats="${WORK_DIR}/bacdive.stats" "${AWK_LIB}"'
    BEGIN {
        while ((getline line < spec) > 0) {
            split(line, f, "\t"); n++
            kind[n] = f[2]; col[n] = f[3] + 0; a[n] = f[4]; b[n] = f[5]; label[n] = f[6]
            if (kind[n] == "class") { m = split(f[4], cl, ","); for (k = 1; k <= m; k++) want[n, cl[k]] = 1 }
        }
        close(spec)
    }
    # classes of the cell "coccus=0.25;rod=0.75" that were requested by filter i, joined with ";"
    function classes_found(cell, i,   parts, m, j, tok, out) {
        if (is_empty(cell)) return ""
        m = split(tolower(cell), parts, ";"); out = ""
        for (j = 1; j <= m; j++) {
            tok = parts[j]; sub(/=.*/, "", tok); gsub(/^[ \t]+|[ \t]+$/, "", tok)
            if ((i, tok) in want) out = out (out == "" ? "" : ";") tok
        }
        return out
    }
    NR == 1 { next }
    {
        if (is_empty($c_taxid)) next
        rank = (is_empty($c_rank) ? "NULL" : $c_rank)
        high = above_species(rank)
        out = $c_taxid OFS (rank == "species" ? "species" : "exact") OFS rank
        for (i = 1; i <= n; i++) {
            v = (kind[i] == "range" ? (any_in_range($(col[i]), a[i], b[i]) ? "1" : "") : classes_found($(col[i]), i))
            if (v != "") { if (high) n_high[i, rank]++; else n_rank[i, rank]++; ranks[rank] = 1 }
            out = out OFS v
        }
        if (!high) print out
    }
    END {
        # per filter, e.g. "BacDive temperature_opt 45:80: 231 NCBI TaxIDs match (200 species, 31 strain)"
        for (i = 1; i <= n; i++) {
            s = ""; t = 0; h = ""
            for (r in ranks) {
                if ((i, r) in n_rank) { s = s (s == "" ? "" : ", ") n_rank[i, r] " " r; t += n_rank[i, r] }
                if ((i, r) in n_high) h = h (h == "" ? "" : ", ") n_high[i, r] " " r
            }
            printf "  %s -> %d NCBI TaxIDs match%s%s\n", label[i], t, (s == "" ? "" : " (" s ")"), \
                (h == "" ? "" : "; above the species rank, not linked to genomes: " h) > stats
        }
    }' "${BACDIVE_TSV}" > "${BD_LINKS_FILE}"
    while IFS= read -r LINE; do
        log "message" "${LINE}"
    done < "${WORK_DIR}/bacdive.stats"
fi

############ GTDB FILTERS ################
#   <key> <TAB> <filter index> <TAB> taxa of the filter found in the GTDB lineage of the genome

GTDB_KEYS_FILE=""
if [[ "${#GTDB_FILTERS[@]}" -gt 0 ]]; then
    GTDB_DIR="${REF_DIR}/GTDB"
    GTDB_TAX_TABLES=()
    for DOMAIN in "${PF_DOMAINS[@]}"; do
        TABLE="$(pf_gtdb_taxonomy_table "${GTDB_DIR}" "${DOMAIN}")"
        [[ -f "${TABLE}" ]] || die "GTDB taxonomy table not found: ${TABLE}. Install the reference databases with: main-prokfinder.sh setup"
        GTDB_TAX_TABLES+=("${TABLE}")
    done
    GTDB_KEYS_FILE="${WORK_DIR}/gtdb_keys.tsv"
    awk -F'\t' -v OFS='\t' -v values="${FILTER_VALUES}" -v spec="${FILTERS_SPEC}" "${AWK_LIB}"'
        BEGIN {
            while ((getline line < spec) > 0) { split(line, f, "\t"); if (f[2] == "GTDB_TAXON") gtdb[f[1]] = 1 }
            close(spec)
            while ((getline line < values) > 0) {
                split(line, f, "\t")
                if (!(f[1] in gtdb)) continue
                v = tolower(f[2]); sub(/^[a-z]__/, "", v)
                want[v] = want[v] SUBSEP f[1] "=" f[2]
            }
            close(values)
        }
        FNR == 1 { delete rank_col; for (i = 1; i <= NF; i++) if ($i ~ /^GTDB_(DOMAIN|PHYLUM|CLASS|ORDER|FAMILY|GENUS|SPECIES)$/) rank_col[i] = 1; next }
        {
            delete found
            for (i in rank_col) {
                v = tolower($i)
                if (!(v in want)) continue
                m = split(substr(want[v], 2), items, SUBSEP)
                for (k = 1; k <= m; k++) { split(items[k], p, "="); found[p[1]] = (p[1] in found ? found[p[1]] ";" : "") substr(items[k], length(p[1]) + 2) }
            }
            for (fi in found) print acc_base($1), fi, found[fi]
        }
    ' "${GTDB_TAX_TABLES[@]}" > "${GTDB_KEYS_FILE}"
fi

############ TAXON NAMES ################
#   <filter index> <TAB> <taxid> <TAB> <name>

TAXON_NAME_FILE=""
if [[ "${HAS_TAXON_NAME}" -eq 1 ]]; then
    TAXON_NAME_FILE="${WORK_DIR}/taxon_name_filters.tsv"
    awk -F'\t' -v OFS='\t' -v spec="${FILTERS_SPEC}" -v names="${WORK_DIR}/taxon_names.tsv" '
        BEGIN {
            while ((getline line < spec) > 0) { split(line, f, "\t"); if (f[2] == "TAXON_NAME") tn[f[1]] = 1 }
            close(spec)
            while ((getline line < names) > 0) { split(line, f, "\t"); ids[tolower(f[2])] = ids[tolower(f[2])] " " f[1] }
            close(names)
        }
        ($1 in tn) { m = split(ids[tolower($2)], t, " "); for (k = 1; k <= m; k++) print $1, t[k], $2 }
    ' "${FILTER_VALUES}" > "${TAXON_NAME_FILE}"
fi

############ NCBI ASSEMBLIES ################
#   Every assembly is checked against every filter; the expression (in postfix notation) decides the selection.

C_ACC="$(require_column "${NCBI_HEADER_FILE}" ASSEMBLY_ACCESSION)"
C_TAXID="$(require_column "${NCBI_HEADER_FILE}" ORGANISM_TAXONOMIC_ID)"
C_ORG="$(require_column "${NCBI_HEADER_FILE}" ORGANISM_NAME)"
C_ISO="$(require_column "${NCBI_HEADER_FILE}" ASSEMBLY_BIOSAMPLE_ISOLATION_SOURCE)"
C_COUNTRY="$(require_column "${NCBI_HEADER_FILE}" ASSEMBLY_BIOSAMPLE_GEOGRAPHIC_LOCATION)"
C_GTYPE="$(require_column "${NCBI_HEADER_FILE}" GENOME_TYPE)"
C_GC="$(require_column "${NCBI_HEADER_FILE}" ASSEMBLY_STATS_GC_PERCENT)"
C_SPECIES_TAXID="$(require_column "${NCBI_HEADER_FILE}" NCBI_SPECIES_TAXID)"
C_STRAIN_TAXID="$(require_column "${NCBI_HEADER_FILE}" NCBI_STRAIN_TAXID)"

# Columns of the NCBI lineage (names and TaxIDs of each rank)
LINEAGE_NAME_COLS="${C_ORG}"
LINEAGE_TAXID_COLS="${C_TAXID}"
for RANK in DOMAIN PHYLUM CLASS ORDER FAMILY GENUS SPECIES STRAIN; do
    LINEAGE_NAME_COLS="${LINEAGE_NAME_COLS},$(require_column "${NCBI_HEADER_FILE}" "NCBI_${RANK}")"
    LINEAGE_TAXID_COLS="${LINEAGE_TAXID_COLS},$(require_column "${NCBI_HEADER_FILE}" "NCBI_${RANK}_TAXID")"
done

SELECTED_TMP="${WORK_DIR}/selected.tsv"

awk -F'\t' -v OFS='\t' \
    -v c_acc="${C_ACC}" -v c_taxid="${C_TAXID}" -v c_iso="${C_ISO}" -v c_country="${C_COUNTRY}" \
    -v c_gtype="${C_GTYPE}" -v c_gc="${C_GC}" -v c_species_taxid="${C_SPECIES_TAXID}" -v c_strain_taxid="${C_STRAIN_TAXID}" \
    -v name_cols="${LINEAGE_NAME_COLS}" -v taxid_cols="${LINEAGE_TAXID_COLS}" \
    -v spec="${FILTERS_SPEC}" -v values="${FILTER_VALUES}" -v bd_file="${BD_LINKS_FILE}" \
    -v gtdb_file="${GTDB_KEYS_FILE}" -v tname_file="${TAXON_NAME_FILE}" \
    -v rpn="${RPN[*]}" -v stats="${WORK_DIR}/selection.stats" "${AWK_LIB}"'
function norm(s) { s = tolower(s); gsub(/[^a-z0-9]+/, " ", s); sub(/ +$/, "", s); return " " s }
function add_unique(list, v) { return (index(";" list ";", ";" v ";") > 0 ? list : list (list == "" ? "" : ";") v) }
BEGIN {
    while ((getline line < spec) > 0) {
        split(line, f, "\t"); n++
        key[n] = f[2]; a1[n] = f[3]; a2[n] = f[4]; sel[n] = f[5]; mname[n] = f[6]; label[n] = f[7]
    }
    close(spec)
    # Values of the filters with a list (nv[i] values; for the text filters, the text to search and the original)
    while ((getline line < values) > 0) {
        split(line, f, "\t"); i = f[1]; v = f[2]
        if (key[i] == "ACCESSION")             set[i, acc_number(v)] = 1
        else if (key[i] == "TAXID")            set[i, v] = v
        else if (key[i] == "NCBI_TAXON")       set[i, tolower(v)] = v
        else if (key[i] == "ISOLATION_SOURCE") { nv[i]++; txt[i, nv[i]] = tolower(v); orig[i, nv[i]] = v }
        else if (key[i] == "COUNTRY")          { nv[i]++; txt[i, nv[i]] = norm(v); orig[i, nv[i]] = v }
    }
    close(values)
    if (tname_file != "") { while ((getline line < tname_file) > 0) { split(line, f, "\t"); set[f[1], f[2]] = (((f[1], f[2]) in set) ? set[f[1], f[2]] ";" : "") f[3] } close(tname_file) }
    if (gtdb_file != "")  { while ((getline line < gtdb_file) > 0)  { split(line, f, "\t"); gtdb[f[2], f[1]] = f[3] } close(gtdb_file) }
    # BacDive: column of each BacDive filter in the links file (4 = first filter)
    nb = 0
    for (i = 1; i <= n; i++) if (key[i] ~ /^BACDIVE_/) bd_col[i] = 3 + (++nb)
    if (bd_file != "") {
        while ((getline line < bd_file) > 0) {
            split(line, f, "\t")
            if (f[2] == "species") bd_species[f[1]] = line; else bd_exact[f[1]] = line
        }
        close(bd_file)
    }
    nr = split(rpn, op, " ")
    n_name_cols = split(name_cols, name_col, ",")
    n_taxid_cols = split(taxid_cols, taxid_col, ",")
}
FNR == 1 {
    if (NR == 1) {
        # each SELECTED_BY column is followed by its column of matched values (filters with several values)
        out = $0
        for (i = 1; i <= n; i++) out = out OFS sel[i] (mname[i] != "" ? OFS mname[i] : "")
        print out
    }
    next
}
{
    n_in++
    key_gtdb = acc_base($c_acc)
    # BacDive row of the assembly: organism/strain TaxID (below species) first, then the species TaxID
    bd = ""
    if (nb > 0) {
        if ($c_taxid in bd_exact) bd = bd_exact[$c_taxid]
        else if ($c_strain_taxid in bd_exact) bd = bd_exact[$c_strain_taxid]
        else if ($c_species_taxid in bd_species) bd = bd_species[$c_species_taxid]
        if (bd != "") split(bd, bdf, "\t")
    }

    for (i = 1; i <= n; i++) {
        m = ""
        k = key[i]
        if (k == "ACCESSION") {
            p = ((i, acc_number($c_acc)) in set)
        } else if (k == "TAXID" || k == "TAXON_NAME") {
            for (j = 1; j <= n_taxid_cols; j++) if ((i, $(taxid_col[j])) in set) m = add_unique(m, set[i, $(taxid_col[j])])
            p = (m != "")
        } else if (k == "NCBI_TAXON") {
            for (j = 1; j <= n_name_cols; j++) if ((i, tolower($(name_col[j]))) in set) m = add_unique(m, set[i, tolower($(name_col[j]))])
            p = (m != "")
        } else if (k == "ISOLATION_SOURCE" || k == "COUNTRY") {
            c = (k == "ISOLATION_SOURCE" ? $c_iso : $c_country)
            if (!is_empty(c)) {
                # isolation source: case-insensitive regular expression; country: start of a word
                if (k == "ISOLATION_SOURCE") {
                    text = tolower(c)
                    for (j = 1; j <= nv[i]; j++) if (text ~ txt[i, j]) m = add_unique(m, orig[i, j])
                } else {
                    text = norm(c)
                    for (j = 1; j <= nv[i]; j++) if (index(text, txt[i, j]) > 0) m = add_unique(m, orig[i, j])
                }
            }
            p = (m != "")
        } else if (k == "GTDB_TAXON") {
            m = (((i, key_gtdb) in gtdb) ? gtdb[i, key_gtdb] : "")
            p = (m != "")
        } else if (k == "GENOME_TYPE") {
            p = (tolower($c_gtype) == a1[i])
        } else if (k == "GC") {
            p = any_in_range($c_gc, a1[i], a2[i])
        } else {
            # BacDive filters
            v = (bd == "" ? "" : bdf[bd_col[i]])
            p = (v != "")
            if (k == "BACDIVE_CLASS") m = v
        }
        pass[i] = p; mval[i] = m
        if (p) n_true[i]++
    }

    # Evaluation of the expression (postfix notation)
    sp = 0
    for (j = 1; j <= nr; j++) {
        if (op[j] == "AND")     { sp--; st[sp] = st[sp] && st[sp + 1] }
        else if (op[j] == "OR") { sp--; st[sp] = st[sp] || st[sp + 1] }
        else                    { st[++sp] = pass[op[j]] }
    }
    if (!st[1]) next
    n_sel++

    out = $0
    for (i = 1; i <= n; i++) out = out OFS (pass[i] ? "TRUE" : "FALSE") (mname[i] != "" ? OFS (mval[i] == "" ? "NULL" : mval[i]) : "")
    print out
}
END {
    printf "NCBI genome assemblies searched\t%d\n", n_in > stats
    for (i = 1; i <= n; i++) printf "  ... matching %s [%s]\t%d\n", label[i], sel[i], n_true[i] > stats
    printf "Selected by the expression\t%d\n", n_sel > stats
}' "${NCBI_TABLES[@]}" > "${SELECTED_TMP}"

while IFS=$'\t' read -r LABEL VALUE; do
    log "message" "${LABEL}: ${VALUE}"
done < "${WORK_DIR}/selection.stats"
TOTAL_SELECTED="$(count_records "${SELECTED_TMP}")"

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

mkdir -p "$(dirname "${OUTPUT_TSV}")"
mv "${SELECTED_TMP}" "${OUTPUT_TSV}"
log "message" "Selected genome metadata (one row per assembly, with the SELECTED_BY_* columns) saved to: ${OUTPUT_TSV}"
