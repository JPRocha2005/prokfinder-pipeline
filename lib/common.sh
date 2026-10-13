#!/usr/bin/env bash
#
# common.sh - functions and paths shared by main-prokfinder.sh and every module.
# This file is sourced, never executed:
#     source "${SCRIPT_DIR}/../lib/common.sh"

# Load only once, even if several files source it
if [[ -n "${PROKFINDER_COMMON_LOADED:-}" ]]; then
    return 0
fi
PROKFINDER_COMMON_LOADED=1

# Associative arrays (declare -A) and "mapfile" need bash 4 or newer
if (( BASH_VERSINFO[0] < 4 )); then
    echo "ERROR: bash 4 or newer is required (found ${BASH_VERSION})." >&2
    exit 1
fi

############ LOGGING ################

# Colors are used only when the messages are shown on a terminal (and NO_COLOR is not set). The decision is made
# once, by the first script (before its output is copied to the log), and inherited by the modules it calls.
if [[ -z "${PF_COLOR:-}" ]]; then
    if [[ -t 2 && -z "${NO_COLOR:-}" ]]; then
        PF_COLOR=1
    else
        PF_COLOR=0
    fi
    export PF_COLOR
fi
if [[ "${PF_COLOR}" -eq 1 ]]; then
    PF_RED=$'\033[1;31m'
    PF_YELLOW=$'\033[1;33m'
    PF_RESET=$'\033[0m'
else
    PF_RED=""
    PF_YELLOW=""
    PF_RESET=""
fi

log() {
    local level="$1"
    local message="$2"
    local timestamp
    timestamp=$(date +"%Y-%m-%d %H:%M:%S")

    case "${level}" in
        "message")
            echo -e "[${timestamp}] [INFO] ${message}"
            ;;
        "warning")
            echo -e "${PF_YELLOW}[${timestamp}] [WARNING] ${message}${PF_RESET}" >&2
            ;;
        "error")
            echo -e "${PF_RED}[${timestamp}] [ERROR] ${message}${PF_RESET}" >&2
            ;;
        *)
            echo -e "[${timestamp}] [LOG] ${message}"
            ;;
    esac
}

# Logs an error and stops the script
die() {
    log "error" "$1"
    exit "${2:-1}"
}

# Temporary directories registered here are removed when the script ends (see _pf_on_exit)
PF_TEMP_DIRS=()
PF_TEE_ACTIVE=0

# Creates a temporary directory and stores its path in the variable named <var_name>.
# (It cannot be used as $(make_temp_dir): the registration would be lost in the subshell.)
# Usage: make_temp_dir <var_name>
make_temp_dir() {
    local DIR
    DIR="$(mktemp -d "${TMPDIR:-/tmp}/prokfinder.XXXXXX")"
    PF_TEMP_DIRS+=("${DIR}")
    printf -v "$1" '%s' "${DIR}"
}

# Single EXIT handler: removes the temporary directories and lets tee finish writing the log
_pf_on_exit() {
    local RC=$?
    local DIR
    for DIR in ${PF_TEMP_DIRS[@]+"${PF_TEMP_DIRS[@]}"}; do
        rm -rf -- "${DIR}"
    done
    if [[ "${PF_TEE_ACTIVE}" -eq 1 ]]; then
        exec 1>&- 2>&-
        wait
    fi
    exit "${RC}"
}
trap _pf_on_exit EXIT

# Copies everything printed by the script (stdout and stderr) to <log_file>, appending to it.
# The terminal keeps the colors; the color codes are removed from the copy written to the log file.
# Usage: start_logging <log_file>
start_logging() {
    mkdir -p "$(dirname "$1")"
    exec 3>&1
    exec > >(tee /dev/fd/3 | sed -u $'s/\033\\[[0-9;]*m//g' >> "$1") 2>&1
    exec 3>&-
    PF_TEE_ACTIVE=1
}

# A module run on its own writes its own log. When it is called by main-prokfinder.sh (PROKFINDER_LOG is set),
# its output already goes to the main log, so nothing is done.
# Usage: setup_module_log <log_file>
setup_module_log() {
    if [[ -z "${PROKFINDER_LOG:-}" ]]; then
        start_logging "$1"
    fi
}

############ FILE AND TABLE HELPERS ################

# Prints the number of data rows of a TSV file (header excluded)
count_records() {
    local LINES
    LINES=$(wc -l < "$1")
    echo $(( LINES > 0 ? LINES - 1 : 0 ))
}

# Prints the number of lines of a file (0 if it does not exist)
count_lines() {
    if [[ -f "$1" ]]; then
        wc -l < "$1"
    else
        echo 0
    fi
}

# Prints the modification time of a file in seconds (0 if it does not exist)
file_mtime() {
    stat -c %Y -- "$1" 2>/dev/null || echo 0
}

# Prints the 1-based position of the first column whose header matches one of the given names.
# The comparison ignores case and every character that is not a letter or a digit, so
# "ASSEMBLY_ACCESSION", "Assembly Accession" and "assemblyaccession" are the same name.
# Prints nothing if there is no match.
# Usage: get_column_index <tsv_file> <name> [alternative_name ...]
get_column_index() {
    local FILE="$1"
    shift
    head -n 1 "${FILE}" | awk -F'\t' -v names="$*" '
    function key(s) { s = tolower(s); gsub(/[^a-z0-9]/, "", s); return s }
    BEGIN { n = split(names, w, " "); for (k = 1; k <= n; k++) wanted[k] = key(w[k]) }
    {
        for (k = 1; k <= n; k++)
            for (i = 1; i <= NF; i++)
                if (key($i) == wanted[k]) { print i; exit }
    }'
}

# Same as get_column_index, but stops the script with a clear message if the column does not exist
# Usage: COLUMN=$(require_column <tsv_file> <name> [alternative_name ...])
require_column() {
    local FILE="$1"
    local INDEX
    shift
    INDEX="$(get_column_index "${FILE}" "$@")"
    if [[ -z "${INDEX}" ]]; then
        log "error" "Column '$1' not found in ${FILE}. Is the file tab-separated?"
        exit 1
    fi
    echo "${INDEX}"
}

# Formats a TSV read from stdin and writes it to stdout:
#   - header: spaces are replaced by "_" (with <upper> = 1 the header is also converted to upper case and every
#     run of characters other than letters/digits becomes "_", e.g. "Assembly BioSample Isolation source" ->
#     ASSEMBLY_BIOSAMPLE_ISOLATION_SOURCE)
#   - cells: trimmed; empty or invalid values ("-", NA, N/A, ND, NaN, none, unknown, missing..., not ...) -> NULL
#   - every row gets the same number of columns as the header (missing cells -> NULL)
# Usage: format_tsv <upper: 0|1> < input.tsv > output.tsv
format_tsv() {
    awk -F'\t' -v OFS='\t' -v upper="${1:-0}" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            h = $i
            gsub(/\r/, "", h)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", h)
            if (upper + 0) { h = toupper(h); gsub(/[^A-Z0-9]+/, "_", h); gsub(/^_+|_+$/, "", h) }
            else gsub(/[[:space:]]+/, "_", h)
            $i = h
        }
        nh = NF
        print
        next
    }
    {
        for (i = 1; i <= nh; i++) {
            v = (i <= NF ? $i : "")
            gsub(/\r/, "", v)
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", v)
            l = tolower(v)
            if (v == "" || l ~ /^-+$/ || l ~ /^(null|na|n\/a|nd|nan|none|unknown|missing)$/ || l ~ /^(not|missing)[ _:-]/)
                v = "NULL"
            $i = v
        }
        NF = nh
        print
    }'
}

# Adds the columns of <table> to <master>. In both files the first column is the assembly accession.
# Columns of <master> with the same name as a column of <table> are replaced, so a stage can be re-run
# and its columns are simply refreshed. Accessions missing from <table> get NULL in the new columns.
# Usage: add_metadata_columns <master_tsv> <table_tsv>
add_metadata_columns() {
    local MASTER="$1"
    local TABLE="$2"

    if [[ ! -s "${TABLE}" ]]; then
        log "warning" "Table ${TABLE} is empty or missing: no column was added to ${MASTER}."
        return 0
    fi

    awk -F'\t' -v OFS='\t' -v table="${TABLE}" '
    BEGIN {
        getline line < table
        nt = split(line, th, "\t")
        for (i = 2; i <= nt; i++) newcol[th[i]] = 1
        while ((getline line < table) > 0) {
            split(line, f, "\t")
            row[f[1]] = line
        }
        close(table)
    }
    NR == 1 {
        nk = 0
        for (i = 1; i <= NF; i++) if (!($i in newcol)) keep[++nk] = i
        out = $(keep[1])
        for (k = 2; k <= nk; k++) out = out OFS $(keep[k])
        for (i = 2; i <= nt; i++) out = out OFS th[i]
        print out
        next
    }
    {
        out = $(keep[1])
        for (k = 2; k <= nk; k++) out = out OFS $(keep[k])
        if ($1 in row) {
            n = split(row[$1], f, "\t")
            for (i = 2; i <= nt; i++) out = out OFS ((i <= n && f[i] != "") ? f[i] : "NULL")
        } else {
            for (i = 2; i <= nt; i++) out = out OFS "NULL"
        }
        print out
    }' "${MASTER}" > "${MASTER}.tmp"
    mv "${MASTER}.tmp" "${MASTER}"
}

# Prints the accessions (first column) of the rows of <tsv> whose column <column_name> has one of the given values
# Usage: select_accessions <tsv> <column_name> <value> [value ...]
select_accessions() {
    local FILE="$1"
    local COLUMN_NAME="$2"
    local COLUMN
    shift 2
    COLUMN="$(require_column "${FILE}" "${COLUMN_NAME}")"
    awk -F'\t' -v c="${COLUMN}" -v values="$*" '
        BEGIN { n = split(values, v, " "); for (i = 1; i <= n; i++) want[v[i]] = 1 }
        NR > 1 && ($c in want) { print $1 }
    ' "${FILE}"
}

# Prints the path of the genome file of <accession> in <dir> (<accession>.fna or <accession>.fna.gz),
# or nothing if there is none
# Usage: genome_file_of <dir> <accession>
genome_file_of() {
    if [[ -f "$1/$2.fna" ]]; then
        echo "$1/$2.fna"
    elif [[ -f "$1/$2.fna.gz" ]]; then
        echo "$1/$2.fna.gz"
    fi
}

# Prints a path as an absolute path (the file or folder must exist)
absolute_path() {
    if [[ -d "$1" ]]; then
        (cd "$1" && pwd)
    else
        echo "$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
    fi
}

# Stops the script if one of the given programs is not in the PATH
# Usage: require_programs <program> [program ...]
require_programs() {
    local PROGRAM
    local MISSING=()
    for PROGRAM in "$@"; do
        if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
            MISSING+=("${PROGRAM}")
        fi
    done
    if [[ "${#MISSING[@]}" -gt 0 ]]; then
        die "Program(s) not found in the PATH: ${MISSING[*]}"
    fi
}

############ PATHS OF THE REFERENCE DATABASE ################
#
# Every module builds the paths of the reference tables with these functions, so the layout of
# Prokfinder/INPUT/REFERENCE-METADATA is defined in a single place.

# Domains handled by Prokfinder (same order everywhere)
# shellcheck disable=SC2034  # used by the modules
PF_DOMAINS=("BACTERIA" "ARCHAEA")

# NCBI (<ncbi_dir> = .../REFERENCE-METADATA/NCBI; <source> = genbank | refseq)
pf_ncbi_genomes_dir()   { echo "$1/$2-NCBI-GENOMES"; }
pf_ncbi_fetch_file()    { echo "$1/$2-NCBI-GENOMES/fetch.txt"; }
pf_ncbi_metadata_dir()  { echo "$1/$2-NCBI-METADATA"; }
pf_ncbi_metadata_table() { echo "$1/$2-NCBI-METADATA/$3-metadata.tsv"; }
pf_ncbi_taxdump_dir()   { echo "$1/TAXDUMP"; }

# BacDive (<bacdive_dir> = .../REFERENCE-METADATA/BACDIVE)
pf_bacdive_brute_dir()  { echo "$1/BACDIVE-DOWNLOAD-BRUTE"; }
pf_bacdive_table()      { echo "$1/BACDIVE-FORMATTED-TSV/BACDIVE_DATA_FORMATTED.tsv"; }

# GTDB (<gtdb_dir> = .../REFERENCE-METADATA/GTDB)
pf_gtdb_prefix() {
    case "$1" in
        BACTERIA) echo "bac120" ;;
        ARCHAEA)  echo "ar53" ;;
        *) die "Unknown domain: $1" ;;
    esac
}
pf_gtdb_taxonomy_table() { echo "$1/$2-GTDB-TAXONOMY/$(pf_gtdb_prefix "$2")_taxonomy_formatted.tsv"; }
pf_gtdb_metadata_table() { echo "$1/$2-GTDB-METADATA/$(pf_gtdb_prefix "$2")_metadata_formatted.tsv"; }
pf_gtdb_version_file()   { echo "$1/VERSION.txt"; }

# Checkpoint written by each download module when it finishes successfully
pf_completed_marker()    { echo "$1/COMPLETED-DOWNLOAD.txt"; }
