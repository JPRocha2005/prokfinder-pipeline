#!/usr/bin/env bash
#
# main-prokfinder.sh - Prokfinder: pipeline for downloading, curating and identifying genomes of cultivated and
# uncultivated prokaryotes.
#
# The command line is organized in sub-commands (run in this order when several are given):
#
#   setup            0. download-reference-db-ncbi.sh / -bacdive.sh / -gtdb.sh   reference databases
#   search-genomes   1. search-genomes-in-reference-db.sh                        genome selection (filters + AND/OR)
#   download-wf      2. download-selected-genomes-from-db.sh                     download (rehydration) of the genomes
#                    3. filter-low-quality-contigs-from-genomes.sh               small contigs removed, sorted, renamed
#                    4. evaluate-and-filter-genome-by-qc.sh                      seqkit + CheckM2 quality filters
#                    5. assign-taxonomy-gtdb.sh                                  GTDB taxonomy (GTDB-Tk optional)
#                    6. build-taxonomic-dir-hierarchy.sh                         Domain/.../Species directory hierarchy
#                    7. dereplicate-genomes-skder.sh                             dereplication inside each taxon (skDER)
#                    8. summarize-final-results.sh                               final tables and graphs
#   run-on-cluster   submit-prokfinder-pbs.sh: the other sub-commands run as a PBS job on a compute node
#
# This script only controls the execution: it reads the command line, checks the programs and the reference
# databases, creates the directory structure and calls the modules (folder modules/). Every step adds its columns to
# a single table, OUTPUT/GENOME-METADATA/genome-metadata.tsv (one row per genome), and can be re-executed: finished
# work is reused.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES_DIR="${SCRIPT_DIR}/modules"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

PROKFINDER_VERSION="0.4.0"
PROG="$(basename "$0")"

############ HELP ############

help_global_options() {
    cat << EOF
Global options (anywhere in the command line):
  -o, --results-dir <dir>   Output directory (default: ./prokfinder-results)
  --refdb-dir <dir>         Reference databases directory (default: \$HOME/prokfinder-refdb)
  -p, --threads <n>         CPUs used by the programs (default: 4; with run-on-cluster, the CPUs of the job)
  --check                   Show which dependencies are installed and exit
  -h, --help                Show the help (of the sub-command, when given after it)
EOF
}

help_main() {
    cat << EOF
Prokfinder ${PROKFINDER_VERSION} - selection, download and curation of prokaryote genomes

Usage: ${PROG} <sub-command> [options]
       ${PROG} --check

Sub-commands (several can be given in the same command; they run in this order):
  setup            Installs the reference databases (NCBI, BacDive and GTDB). Run it once, before the others
  search-genomes   Selects the genomes that match the selection filters
  download-wf      Downloads the selected genomes, filters them by quality, assigns the GTDB taxonomy and removes
                   the redundant genomes
  run-on-cluster   Runs the other sub-commands as a job on a cluster (PBS)

$(help_global_options)

Examples:
  ${PROG} setup
  ${PROG} search-genomes -T 45:80 AND -O rod,coccus
  ${PROG} -o thermophiles search-genomes -T 45:80 AND -O rod,coccus download-wf --min-completeness 90
  ${PROG} run-on-cluster search-genomes --temp-growth 45:80 AND -i "hot spring"

Help of a sub-command: ${PROG} <sub-command> -h
EOF
}

help_setup() {
    cat << EOF
Usage: ${PROG} setup [options]

Downloads and installs the reference databases (NCBI, BacDive and GTDB) in the reference databases directory.
It takes some hours and is made only once: the databases are shared by every run.

Options:
  --force-refdb              Download and install the databases again (the current ones are kept until the new
                             ones are complete)
  --ncbi-complete-metadata   Also write a table with every NCBI metadata field (slow; not used by the search)

$(help_global_options)
EOF
}

help_search() {
    cat << EOF
Usage: ${PROG} search-genomes <filter> [AND|OR <filter> ...] [options]

Selects, in the NCBI tables, the genomes that match the selection filters and writes
OUTPUT/GENOME-METADATA/genome-metadata.tsv in the output directory.

How to write the filters:
  - Each flag below is a filter. Filters are joined with AND or OR:   -T 45:80 AND -c Brazil
  - Several values of one flag, separated by commas, match ANY value: -c Brazil,Iceland
  - Values with spaces go between quotes:                             -i "hot spring","hydrothermal vent"
  - Quoted parentheses group filters:                                 -T 45:80 AND "(" -O rod OR -M anaerobe ")"
  - Ranges are <min>:<max>; one side may be empty:                    -T 45:80, -T 70:
  - List files (-f, -A, -N, -I, -C) have one value per line

NCBI filters:
  -t, --taxid <id>                     NCBI TaxID of any rank (e.g., 1239 selects every Bacillota genome)
  -f, --taxid-file <file>              List of NCBI TaxIDs
  -a, --accession <acc>                Accession with version, GCA_ or GCF_ (e.g., GCA_000005845.2)
  -A, --accession-file <file>          List of accessions
  --ncbi-taxon <name>                  NCBI taxon name, any rank (e.g., Thermus)
  -N, --ncbi-taxon-file <file>         List of NCBI taxon names (e.g., Thermus aquaticus)
  -i, --isolation-source <text>        Text in the isolation source (e.g., "hot spring")
  -I, --isolation-source-file <file>   List of isolation sources
  -c, --country <name>                 Country of the geographic location (e.g., Brazil)
  -C, --country-file <file>            List of countries
  -g, --genome-type <type>             isolate_genome or mag
  -G, --gc <min:max>                   GC content in % (e.g., 40:65)

BacDive filters (phenotype of the species or strain):
  -T, --temp-opt <min:max>             Optimal growth temperature in °C (e.g., 45:80)
  --temp-growth <min:max>              Growth temperature in °C (known for many more species than -T)
  -P, --ph-opt <min:max>               Optimal pH (e.g., 5.5:8.5)
  --ph-growth <min:max>                Growth pH
  -S, --salt-opt <min:max>             Optimal NaCl concentration in % (e.g., 0:5)
  --salt-growth <min:max>              Growth NaCl concentration in %
  -O, --cell-shape <shape>             rod, coccus, ovoid, curved, spiral, filament, pleomorphic or other
  -M, --oxygen <tolerance>             aerobe, anaerobe, facultative or microaerophile
  --gram-stain <result>                positive or negative
  --motility <yes|no>                  Motile cells
  --spore-formation <yes|no>           Forms spores
  --energy-source <source>             chemotroph or phototroph
  --electron-donor <donor>             lithotroph or organotroph
  --carbon-source <source>             autotroph or heterotroph

GTDB filters:
  --gtdb-taxon <name>                  GTDB taxon name, any rank (e.g., g__Thermus or Thermus)

Options:
  -s, --assembly-source <src>          NCBI table searched: genbank (default) or refseq
  --force-selection                    Search again, even if the filters did not change

$(help_global_options)

Examples:
  ${PROG} search-genomes -T 45:80 AND -O rod,coccus AND -g isolate_genome
  ${PROG} search-genomes -N taxa.txt OR "(" --temp-growth 70: AND -i "hot spring","hydrothermal vent" ")"
  ${PROG} search-genomes --gtdb-taxon g__Thermus AND -c Brazil,Iceland -s refseq
EOF
}

help_download() {
    cat << EOF
Usage: ${PROG} download-wf [options]

Downloads the genomes selected by search-genomes, removes the small contigs, filters the genomes by quality,
assigns the GTDB taxonomy and removes the redundant genomes of each taxon. Writes the final genomes, tables and
graphs in the output directory (the same one used by search-genomes).

Download:
  -j, --download-threads <n>       Concurrent genome downloads, 1-30 (default: 10)
  -z, --gzip                       Keep the genome files compressed (.fna.gz)

Contig filtering:
  --min-contig-length <bp>         Contigs shorter than this are removed (default: 1000)

Genome quality (seqkit and CheckM2):
  --max-contigs <n>                Maximum number of contigs (default: 1000)
  --max-n-bases <n>                Maximum number of ambiguous bases, N (default: 100000)
  --min-n50 <bp>                   Minimum N50 (default: 5000)
  --min-completeness <pct>         Minimum completeness (default: 75)
  --max-contamination <pct>        Maximum contamination (default: 5)

Taxonomy:
  --gtdbtk                         Classify with GTDB-Tk the genomes that are not in the GTDB tables
                                   (needs more than 200 GB of RAM)

Dereplication (skDER, inside each taxon):
  --derep-ani <pct>                ANI above which two genomes are redundant (default: 99.5)
  --derep-mode <mode>              greedy (default) or dynamic
  --filter-mge                     Mask the mobile genetic elements (PhiSpy) before comparing the genomes

Re-running:
  --force-curation                 Redo the curation from the downloaded genomes (by default, finished work is reused)

$(help_global_options)

Examples:
  ${PROG} download-wf --min-completeness 90 --max-contamination 2
  ${PROG} search-genomes -T 45:80 AND -O rod download-wf --derep-ani 99 --filter-mge -z
EOF
}

help_cluster() {
    cat << EOF
Usage: ${PROG} run-on-cluster [options] <sub-commands>

Runs the other sub-commands as a job on a cluster node, with the PBS queue manager. The job script and its log are
written to PBS-JOBS/ in the output directory.

Options:
  --queue <queue>         Submission queue (default: roteamento):
                            qtime       tests and short runs (walltime of at most 72:00:00)
                            roteamento  common runs (at most 24024 MB of RAM per job)
                            high_mem    runs that need a lot of RAM (25 GB or more; e.g., download-wf --gtdbtk)
  --mem <ram>             RAM of the job: 20gb, 20000mb or only the number in GB (default: 20gb)
  --cpus <n>              CPUs of the job (default: 20; give setup at least 4)
  --walltime <HH:MM:SS>   Maximum run time of the job (default: 2160:00:00; 72:00:00 with qtime)
  --test                  Submit only a short job that shows the version of every program
  --dry-run               Only write the job script (do not submit it)

$(help_global_options)

Examples:
  ${PROG} run-on-cluster setup
  ${PROG} run-on-cluster --queue qtime --walltime 24:00:00 search-genomes -T 45:80 AND -O rod
  ${PROG} run-on-cluster --queue high_mem --mem 256gb --cpus 32 download-wf --gtdbtk
  ${PROG} run-on-cluster --test
EOF
}

show_help() {
    case "$1" in
        setup)          help_setup ;;
        search-genomes) help_search ;;
        download-wf)    help_download ;;
        run-on-cluster) help_cluster ;;
        *)              help_main ;;
    esac
    exit 0
}

# Command line errors: only the error and a hint are printed (the help is shown only with -h/--help)
cli_error() {
    log "error" "$1"
    echo "Use ${PROG} ${2:+$2 }-h (or --help) to see the options." >&2
    exit 1
}

############ READING THE COMMAND LINE ############

SUBCOMMANDS=" setup search-genomes download-wf run-on-cluster "

# Options of each section: flags (no value) and options with a value
GLOBAL_FLAGS=" -h --help --check "
GLOBAL_VALUES=" -o --results-dir --refdb-dir -p --threads "
SETUP_FLAGS=" --force-refdb --ncbi-complete-metadata "
SETUP_VALUES=" "
SEARCH_FLAGS=" --force-selection "
SEARCH_VALUES=" -s --assembly-source "
SEARCH_FILTERS=" -t --taxid -f --taxid-file -a --accession -A --accession-file --ncbi-taxon -N --ncbi-taxon-file
 -i --isolation-source -I --isolation-source-file -c --country -C --country-file -g --genome-type -G --gc
 -T --temp-opt -P --ph-opt -S --salt-opt --temp-growth --ph-growth --salt-growth -O --cell-shape -M --oxygen
 --gram-stain --motility --spore-formation --energy-source --electron-donor --carbon-source --gtdb-taxon "
SEARCH_FILTERS="${SEARCH_FILTERS//$'\n'/ }"
DOWNLOAD_FLAGS=" -z --gzip --gtdbtk --filter-mge --force-curation "
DOWNLOAD_VALUES=" -j --download-threads --min-contig-length --max-contigs --max-n-bases --min-n50 --min-completeness
 --max-contamination --derep-mode --derep-ani "
DOWNLOAD_VALUES="${DOWNLOAD_VALUES//$'\n'/ }"
CLUSTER_FLAGS=" --test --dry-run "
CLUSTER_VALUES=" --queue --mem --cpus --walltime "

in_list() { [[ "$2" == *" $1 "* ]]; }

# Prints the sub-command that owns an option ("global" or nothing if unknown)
option_section() {
    if in_list "$1" "${GLOBAL_FLAGS}${GLOBAL_VALUES}"; then echo "global"
    elif in_list "$1" "${SETUP_FLAGS}${SETUP_VALUES}"; then echo "setup"
    elif in_list "$1" "${SEARCH_FLAGS}${SEARCH_VALUES}${SEARCH_FILTERS}"; then echo "search-genomes"
    elif in_list "$1" "${DOWNLOAD_FLAGS}${DOWNLOAD_VALUES}"; then echo "download-wf"
    elif in_list "$1" "${CLUSTER_FLAGS}${CLUSTER_VALUES}"; then echo "run-on-cluster"
    fi
}
takes_value() { in_list "$1" "${GLOBAL_VALUES}${SETUP_VALUES}${SEARCH_VALUES}${SEARCH_FILTERS}${DOWNLOAD_VALUES}${CLUSTER_VALUES}"; }
is_known_option() { [[ -n "$(option_section "$1")" ]]; }

# Defaults
PROKFINDER_DIR="prokfinder-results"
REFDB_DIR="${HOME}/prokfinder-refdb"
THREADS=4
THREADS_GIVEN=0
CHECK_DEPS=0
FORCE_REFDB=0
NCBI_COMPLETE=0
FORCE_SELECTION=0
ASSEMBLY_SOURCE="genbank"
DOWNLOAD_THREADS=10
GZIP=0
MIN_CONTIG_LENGTH=1000
QC_ARGS=()
USE_GTDBTK=0
DEREP_MODE="greedy"
DEREP_ANI="99.5"
FILTER_MGE=0
FORCE_CURATION=0
PBS_ARGS=()
declare -A RUN=()

# Selection filters: one id per filter (F1, F2...), its lines in the query file and the expression tokens
QUERY_LINES=()
EXPR_TOKENS=()
N_FILTERS=0
EXPR_PREV="start"     # start | filter | operator | open | close
EXPR_PREV_FLAG=""
EXPR_DEPTH=0

# Adds one value to the current filter
add_value() {
    QUERY_LINES+=("FILTER"$'\t'"F${N_FILTERS}"$'\t'"$1"$'\t'"$2")
}
# Adds one value per comma-separated item
add_values_csv() {
    local ITEM
    local -a ITEMS
    IFS=',' read -ra ITEMS <<< "$2"
    for ITEM in "${ITEMS[@]}"; do
        ITEM="$(echo "${ITEM}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
        [[ -z "${ITEM}" ]] || add_value "$1" "${ITEM}"
    done
}
# Adds one value per line of a list file (every \r removed, lines trimmed, empty lines skipped)
add_values_file() {
    local LINE
    local N=0
    [[ -f "$2" ]] || cli_error "List file not found: '$2'" search-genomes
    while IFS= read -r LINE || [[ -n "${LINE}" ]]; do
        LINE="$(echo "${LINE//$'\r'/}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
        if [[ -n "${LINE}" ]]; then
            add_value "$1" "${LINE}"
            N=$((N + 1))
        fi
    done < "$2"
    [[ "${N}" -gt 0 ]] || cli_error "No valid entry found in the list file '$2'." search-genomes
}

# A new filter: AND/OR must separate it from the previous filter
start_filter() {
    if [[ "${EXPR_PREV}" == "filter" || "${EXPR_PREV}" == "close" ]]; then
        cli_error "AND or OR is missing between ${EXPR_PREV_FLAG} and $1. Two filters must be joined with AND or OR (e.g., -T 45:80 AND -O rod)." search-genomes
    fi
    N_FILTERS=$((N_FILTERS + 1))
    EXPR_TOKENS+=("F${N_FILTERS}")
    EXPR_PREV="filter"
    EXPR_PREV_FLAG="$1"
}

# Expression tokens: AND, OR, ( and )
expression_token() {
    local TOKEN="$1"
    case "${TOKEN}" in
        AND|OR)
            if [[ "${EXPR_PREV}" != "filter" && "${EXPR_PREV}" != "close" ]]; then
                cli_error "${TOKEN} must be between two selection filters (e.g., -T 45:80 ${TOKEN} -O rod)." search-genomes
            fi
            EXPR_PREV="operator"
            ;;
        "(")
            if [[ "${EXPR_PREV}" == "filter" || "${EXPR_PREV}" == "close" ]]; then
                cli_error "AND or OR is missing between ${EXPR_PREV_FLAG} and '('." search-genomes
            fi
            EXPR_DEPTH=$((EXPR_DEPTH + 1))
            EXPR_PREV="open"
            ;;
        ")")
            if [[ "${EXPR_PREV}" != "filter" && "${EXPR_PREV}" != "close" ]]; then
                cli_error "')' must close a group of filters (e.g., \"(\" -O rod OR -M anaerobe \")\")." search-genomes
            fi
            EXPR_DEPTH=$((EXPR_DEPTH - 1))
            [[ "${EXPR_DEPTH}" -ge 0 ]] || cli_error "')' without '('." search-genomes
            EXPR_PREV="close"
            EXPR_PREV_FLAG="')'"
            ;;
    esac
    EXPR_TOKENS+=("${TOKEN}")
}

COMMAND_LINE="$0 $*"
SECTION=""
JOB_ARGS=()   # the command line without the run-on-cluster sub-command and its options (for the PBS job)

while [[ $# -gt 0 ]]; do
    ARG="$1"

    # "--option=value" is accepted as "--option value"
    if [[ "${ARG}" == --*=* ]]; then
        shift
        set -- "${ARG%%=*}" "${ARG#*=}" "$@"
        ARG="$1"
    fi

    # Sub-commands
    if in_list "${ARG}" "${SUBCOMMANDS}"; then
        [[ -z "${RUN[${ARG}]:-}" ]] || cli_error "The sub-command ${ARG} was given more than once."
        RUN["${ARG}"]=1
        SECTION="${ARG}"
        [[ "${ARG}" == "run-on-cluster" ]] || JOB_ARGS+=("${ARG}")
        shift
        continue
    fi

    # AND, OR and parentheses (only between the filters of search-genomes)
    if [[ "${ARG^^}" == "AND" || "${ARG^^}" == "OR" || "${ARG}" == "(" || "${ARG}" == ")" ]]; then
        [[ "${SECTION}" == "search-genomes" ]] || cli_error "'${ARG}' is only used between the selection filters of search-genomes." "${SECTION}"
        TOKEN="${ARG}"
        [[ "${ARG}" == "(" || "${ARG}" == ")" ]] || TOKEN="${ARG^^}"
        expression_token "${TOKEN}"
        JOB_ARGS+=("${TOKEN}")
        shift
        continue
    fi

    if [[ "${ARG}" == "-h" || "${ARG}" == "--help" ]]; then
        show_help "${SECTION}"
    fi

    OWNER="$(option_section "${ARG}")"
    if [[ -z "${OWNER}" ]]; then
        if [[ "${ARG}" == -* ]]; then
            cli_error "Unknown option: '${ARG}'." "${SECTION}"
        fi
        if [[ -z "${SECTION}" ]]; then
            cli_error "Unknown sub-command: '${ARG}'. Sub-commands: setup, search-genomes, download-wf, run-on-cluster."
        fi
        cli_error "Unexpected argument '${ARG}': it is neither an option nor the value of one. Several values of the same option are separated by commas, without spaces (e.g., -t 562,1280), each value with spaces between its own double quotes (e.g., -i \"hot spring\",geyser), and two filters are joined with AND or OR." "${SECTION}"
    fi
    if [[ "${OWNER}" != "global" && "${OWNER}" != "${SECTION}" ]]; then
        if [[ -z "${SECTION}" ]]; then
            cli_error "Option ${ARG} belongs to the sub-command ${OWNER}: write it after the sub-command (e.g., ${PROG} ${OWNER} ${ARG} ...)." "${OWNER}"
        fi
        cli_error "Option ${ARG} belongs to the sub-command ${OWNER}, not to ${SECTION}." "${OWNER}"
    fi

    VAL=""
    if takes_value "${ARG}"; then
        # The value is missing when the option is the last argument or is followed by another option or sub-command
        # (a negative number such as "-T -5:10" is a value, not an option)
        if [[ $# -lt 2 ]] || is_known_option "$2" || in_list "$2" "${SUBCOMMANDS}"; then
            cli_error "Option ${ARG} requires a value." "${OWNER/global/}"
        fi
        VAL="$2"
        shift 2
    else
        shift
    fi

    # Arguments kept for the PBS job (everything except the options of run-on-cluster)
    if [[ "${OWNER}" == "run-on-cluster" ]]; then
        PBS_ARGS+=("${ARG}" ${VAL:+"${VAL}"})
    else
        JOB_ARGS+=("${ARG}")
        takes_value "${ARG}" && JOB_ARGS+=("${VAL}")
    fi

    case "${ARG}" in
        # Global
        -o|--results-dir)             PROKFINDER_DIR="${VAL%/}" ;;
        --refdb-dir)                  REFDB_DIR="${VAL%/}" ;;
        -p|--threads)                 THREADS="${VAL}"; THREADS_GIVEN=1 ;;
        --check)                      CHECK_DEPS=1 ;;
        # setup
        --force-refdb)                FORCE_REFDB=1 ;;
        --ncbi-complete-metadata)     NCBI_COMPLETE=1 ;;
        # search-genomes: options
        -s|--assembly-source)         ASSEMBLY_SOURCE="${VAL,,}"; QUERY_LINES+=("ASSEMBLY_SOURCE"$'\t'"${VAL,,}") ;;
        --force-selection)            FORCE_SELECTION=1 ;;
        # search-genomes: filters
        -t|--taxid)                   start_filter "${ARG}"; add_values_csv TAXID "${VAL}" ;;
        -f|--taxid-file)              start_filter "${ARG}"; add_values_file TAXID "${VAL}" ;;
        -a|--accession)               start_filter "${ARG}"; add_values_csv ACCESSION "${VAL}" ;;
        -A|--accession-file)          start_filter "${ARG}"; add_values_file ACCESSION "${VAL}" ;;
        --ncbi-taxon)                 start_filter "${ARG}"; add_values_csv NCBI_TAXON "${VAL}" ;;
        -N|--ncbi-taxon-file)         start_filter "${ARG}"; add_values_file TAXON_NAME "${VAL}" ;;
        -i|--isolation-source)        start_filter "${ARG}"; add_values_csv ISOLATION_SOURCE "${VAL}" ;;
        -I|--isolation-source-file)   start_filter "${ARG}"; add_values_file ISOLATION_SOURCE "${VAL}" ;;
        -c|--country)                 start_filter "${ARG}"; add_values_csv COUNTRY "${VAL}" ;;
        -C|--country-file)            start_filter "${ARG}"; add_values_file COUNTRY "${VAL}" ;;
        -g|--genome-type)             start_filter "${ARG}"; add_value GENOME_TYPE "${VAL}" ;;
        -G|--gc)                      start_filter "${ARG}"; add_value GC "${VAL}" ;;
        -T|--temp-opt)                start_filter "${ARG}"; add_value BACDIVE_RANGE "temperature_opt=${VAL}" ;;
        -P|--ph-opt)                  start_filter "${ARG}"; add_value BACDIVE_RANGE "ph_opt=${VAL}" ;;
        -S|--salt-opt)                start_filter "${ARG}"; add_value BACDIVE_RANGE "halophily_opt=${VAL}" ;;
        --temp-growth)                start_filter "${ARG}"; add_value BACDIVE_RANGE "temperature_preference=${VAL}" ;;
        --ph-growth)                  start_filter "${ARG}"; add_value BACDIVE_RANGE "ph_preference=${VAL}" ;;
        --salt-growth)                start_filter "${ARG}"; add_value BACDIVE_RANGE "halophily_preference=${VAL}" ;;
        -O|--cell-shape)              start_filter "${ARG}"; add_value BACDIVE_CLASS "cell_shape=${VAL}" ;;
        -M|--oxygen)                  start_filter "${ARG}"; add_value BACDIVE_CLASS "oxygen_tolerance=${VAL}" ;;
        --gram-stain)                 start_filter "${ARG}"; add_value BACDIVE_CLASS "gram_stain=${VAL}" ;;
        --motility)                   start_filter "${ARG}"; add_value BACDIVE_CLASS "motility=${VAL}" ;;
        --spore-formation)            start_filter "${ARG}"; add_value BACDIVE_CLASS "spore_formation=${VAL}" ;;
        --energy-source)              start_filter "${ARG}"; add_value BACDIVE_CLASS "energy_source=${VAL}" ;;
        --electron-donor)             start_filter "${ARG}"; add_value BACDIVE_CLASS "electron_donor=${VAL}" ;;
        --carbon-source)              start_filter "${ARG}"; add_value BACDIVE_CLASS "carbon_source=${VAL}" ;;
        --gtdb-taxon)                 start_filter "${ARG}"; add_values_csv GTDB_TAXON "${VAL}" ;;
        # download-wf
        -j|--download-threads)        DOWNLOAD_THREADS="${VAL}" ;;
        -z|--gzip)                    GZIP=1 ;;
        --min-contig-length)          MIN_CONTIG_LENGTH="${VAL}" ;;
        --max-contigs|--max-n-bases|--min-n50|--min-completeness|--max-contamination) QC_ARGS+=("${ARG}" "${VAL}") ;;
        --gtdbtk)                     USE_GTDBTK=1 ;;
        --derep-mode)                 DEREP_MODE="${VAL,,}" ;;
        --derep-ani)                  DEREP_ANI="${VAL}" ;;
        --filter-mge)                 FILTER_MGE=1 ;;
        --force-curation)             FORCE_CURATION=1 ;;
    esac
done

# Reference databases (shared by every run and every results folder)
REF_DIR="${REFDB_DIR}"
NCBI_DIR="${REF_DIR}/NCBI"
BACDIVE_DIR="${REF_DIR}/BACDIVE"
GTDB_DIR="${REF_DIR}/GTDB"
SEARCH_CHECKPOINT="${PROKFINDER_DIR}/SEARCH-GENOMES-COMPLETED.txt"

# State of the reference databases: REFDB_MISSING[<DB>] = 1 (not installed) or REFDB_OUTDATED[<DB>] = 1 (installed by
# an older Prokfinder, whose tables have another format)
declare -A REFDB_MISSING=()
declare -A REFDB_OUTDATED=()
check_refdb_state() {
    local DB DB_DIR FORMAT_VAR
    REFDB_MISSING=()
    REFDB_OUTDATED=()
    for DB in NCBI BACDIVE GTDB; do
        case "${DB}" in
            NCBI)    DB_DIR="${NCBI_DIR}" ;;
            BACDIVE) DB_DIR="${BACDIVE_DIR}" ;;
            GTDB)    DB_DIR="${GTDB_DIR}" ;;
        esac
        FORMAT_VAR="PF_REFDB_FORMAT_${DB}"
        if [[ ! -f "$(pf_completed_marker "${DB_DIR}")" ]]; then
            REFDB_MISSING["${DB}"]=1
        elif ! pf_refdb_format_ok "${DB_DIR}" "${!FORMAT_VAR}"; then
            REFDB_OUTDATED["${DB}"]=1
        fi
    done
}
check_refdb_state

############ DEPENDENCIES ############
#   setup and search-genomes (the genome search) only need the NCBI Datasets programs and taxonkit, besides common
#   Unix tools; download-wf (the rest of the pipeline) also needs the genome analysis programs. The lists below are
#   used both by the checks made before a run and by --check.

# <program>|<how to install it>
DEPS_SEARCH=(
    "datasets|conda install -c conda-forge ncbi-datasets-cli"
    "dataformat|conda install -c conda-forge ncbi-datasets-cli (installed together with datasets)"
    "taxonkit|conda install -c bioconda taxonkit"
)
DEPS_DOWNLOAD=(
    "seqkit|conda install -c bioconda seqkit"
    "checkm2|conda install -c bioconda checkm2"
    "skder|conda install -c bioconda skder"
)
DEPS_OPTIONAL=(
    "gtdbtk|conda install -c bioconda gtdbtk (only for download-wf --gtdbtk)"
    "Rscript|conda install -c conda-forge r-base (graphs of the final results)"
)
# Unix tools used by setup and the other steps (usually already installed)
DEPS_SYSTEM=(awk sort gzip tar curl unzip md5sum jq gawk)

install_hint() {
    local ENTRY
    for ENTRY in "${DEPS_SEARCH[@]}" "${DEPS_DOWNLOAD[@]}" "${DEPS_OPTIONAL[@]}"; do
        if [[ "${ENTRY%%|*}" == "$1" ]]; then
            echo "${ENTRY#*|}"
            return 0
        fi
    done
    echo "install it with the package manager of the system (e.g., apt install $1) or with conda"
}

# Prints the first line printed by a program to get its version
tool_version() {
    local VERSION
    case "$1" in
        seqkit|dataformat|taxonkit) VERSION="$("$1" version 2>&1 | head -n 1)" || true ;;
        *)                          VERSION="$("$1" --version 2>&1 | head -n 1)" || true ;;
    esac
    [[ "${VERSION}" != "undefined" ]] || VERSION=""
    echo "${VERSION:-(version not reported)}"
}

# Prints what is wrong with CHECKM2DB / GTDBTK_DATA_PATH (nothing when they are fine)
checkm2db_problem() {
    if [[ -z "${CHECKM2DB:-}" ]]; then
        echo "CHECKM2DB is not set. It must be the path of the CheckM2 database file, e.g. export CHECKM2DB=/path/CheckM2_database/uniref100.KO.1.dmnd (download it with: checkm2 database --download; on the cluster, the checkm2 module sets it)."
    elif [[ ! -f "${CHECKM2DB}" || ! -r "${CHECKM2DB}" ]]; then
        echo "CHECKM2DB must point to the database file itself (uniref100.KO.1.dmnd); not found or not readable: ${CHECKM2DB}"
    fi
}
gtdbtk_data_problem() {
    if [[ -z "${GTDBTK_DATA_PATH:-}" || ! -d "${GTDBTK_DATA_PATH}" ]]; then
        echo "GTDBTK_DATA_PATH is not set or is not a folder. It must point to the GTDB-Tk reference package (on the cluster, the gtdbtk module sets it)."
    fi
}

# --check: lists the dependencies of the genome search and of the complete pipeline, and the reference databases
show_dependencies() {
    local GREEN="" RED="" YELLOW="" RESET="" ENTRY NAME TOOL DB
    local MISSING_SEARCH=0 MISSING_FULL=0 REFDB_NOTE=""
    local -a MISSING_SYSTEM=()
    if [[ -t 1 && -z "${NO_COLOR:-}" ]]; then
        GREEN=$'\033[1;32m'; RED=$'\033[1;31m'; YELLOW=$'\033[1;33m'; RESET=$'\033[0m'
    fi
    dep_ok()      { printf "  ${GREEN}%-9s${RESET} %-17s %s\n" "[OK]" "$1" "$2"; }
    dep_missing() { printf "  ${RED}%-9s${RESET} %-17s %s\n" "[MISSING]" "$1" "$2"; }
    dep_absent()  { printf "  ${YELLOW}%-9s${RESET} %-17s %s\n" "[--]" "$1" "$2"; }
    # Usage: program_line <program|hint> <dep_missing|dep_absent>; returns 1 when the program is not found
    program_line() {
        NAME="${1%%|*}"
        if command -v "${NAME}" >/dev/null 2>&1; then
            dep_ok "${NAME}" "$(tool_version "${NAME}")"
        else
            "$2" "${NAME}" "install: ${1#*|}"
            return 1
        fi
    }

    echo "Prokfinder ${PROKFINDER_VERSION} - dependencies"
    echo
    echo "Genome search (setup and search-genomes):"
    for ENTRY in "${DEPS_SEARCH[@]}"; do
        program_line "${ENTRY}" dep_missing || MISSING_SEARCH=$((MISSING_SEARCH + 1))
    done
    for TOOL in "${DEPS_SYSTEM[@]}"; do
        command -v "${TOOL}" >/dev/null 2>&1 || MISSING_SYSTEM+=("${TOOL}")
    done
    if [[ "${#MISSING_SYSTEM[@]}" -eq 0 ]]; then
        dep_ok "Unix tools" "${DEPS_SYSTEM[*]}"
    else
        dep_missing "Unix tools" "${MISSING_SYSTEM[*]} (install with the package manager of the system or with conda)"
        MISSING_SEARCH=$((MISSING_SEARCH + ${#MISSING_SYSTEM[@]}))
    fi

    echo
    echo "Complete pipeline (download-wf), in addition:"
    for ENTRY in "${DEPS_DOWNLOAD[@]}"; do
        program_line "${ENTRY}" dep_missing || MISSING_FULL=$((MISSING_FULL + 1))
    done
    if [[ -z "$(checkm2db_problem)" ]]; then
        dep_ok "CHECKM2DB" "${CHECKM2DB}"
    else
        dep_missing "CHECKM2DB" "not set or not a file: path of the CheckM2 database (uniref100.KO.1.dmnd; download it with: checkm2 database --download)"
        MISSING_FULL=$((MISSING_FULL + 1))
    fi

    echo
    echo "Optional:"
    for ENTRY in "${DEPS_OPTIONAL[@]}"; do
        program_line "${ENTRY}" dep_absent || true
    done
    if [[ -z "$(gtdbtk_data_problem)" ]]; then
        dep_ok "GTDBTK_DATA_PATH" "${GTDBTK_DATA_PATH}"
    else
        dep_absent "GTDBTK_DATA_PATH" "not set (only for download-wf --gtdbtk)"
    fi
    if [[ -n "${NCBI_API_KEY:-}" ]]; then
        dep_ok "NCBI_API_KEY" "set (datasets: up to 10 requests per second)"
    else
        dep_absent "NCBI_API_KEY" "not set (datasets: up to 5 requests per second)"
    fi

    echo
    echo "Reference databases (${REF_DIR}):"
    for DB in NCBI BACDIVE GTDB; do
        if [[ -n "${REFDB_MISSING[${DB}]:-}" ]]; then
            dep_absent "${DB}" "not installed"
            REFDB_NOTE=" (install the reference databases first: ${PROG} setup)"
        elif [[ -n "${REFDB_OUTDATED[${DB}]:-}" ]]; then
            dep_absent "${DB}" "installed by an older Prokfinder"
            REFDB_NOTE=" (update the reference databases first: ${PROG} setup)"
        else
            dep_ok "${DB}" "installed"
        fi
    done

    echo
    if [[ "${MISSING_SEARCH}" -eq 0 ]]; then
        echo "Genome search (setup, search-genomes): ${GREEN}ready${RESET}${REFDB_NOTE}"
    else
        echo "Genome search (setup, search-genomes): ${RED}${MISSING_SEARCH} missing${RESET}"
    fi
    if [[ $((MISSING_SEARCH + MISSING_FULL)) -eq 0 ]]; then
        echo "Complete pipeline (download-wf): ${GREEN}ready${RESET}"
    else
        echo "Complete pipeline (download-wf): ${RED}$((MISSING_SEARCH + MISSING_FULL)) missing${RESET}"
    fi
    if command -v qsub >/dev/null 2>&1 && ! command -v seqkit >/dev/null 2>&1; then
        echo "On the cluster, the programs are loaded by the environment modules of the job: use run-on-cluster --test."
    fi
    [[ $((MISSING_SEARCH + MISSING_FULL)) -eq 0 ]]
}

if [[ "${CHECK_DEPS}" -eq 1 ]]; then
    STATUS=0
    show_dependencies || STATUS=1
    exit "${STATUS}"
fi

############ CHECKING THE COMMAND LINE ############

if [[ "${#RUN[@]}" -eq 0 ]]; then
    cli_error "No sub-command given. Sub-commands: setup, search-genomes, download-wf, run-on-cluster."
fi
if [[ -n "${RUN[run-on-cluster]:-}" && "${#RUN[@]}" -eq 1 && "${PBS_ARGS[*]:-}" != *--test* ]]; then
    cli_error "run-on-cluster needs another sub-command to run on the cluster (setup, search-genomes and/or download-wf), or --test." run-on-cluster
fi

if [[ -n "${RUN[search-genomes]:-}" ]]; then
    [[ "${N_FILTERS}" -gt 0 ]] || cli_error "search-genomes needs at least one selection filter." search-genomes
    if [[ "${EXPR_PREV}" == "operator" || "${EXPR_PREV}" == "open" ]]; then
        cli_error "The selection filters cannot end with AND, OR or '('." search-genomes
    fi
    [[ "${EXPR_DEPTH}" -eq 0 ]] || cli_error "'(' without ')' in the selection filters." search-genomes
    [[ "${ASSEMBLY_SOURCE}" == "genbank" || "${ASSEMBLY_SOURCE}" == "refseq" ]] \
        || cli_error "Invalid -s/--assembly-source '${ASSEMBLY_SOURCE}'. Must be genbank or refseq." search-genomes
fi

[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || cli_error "Invalid number of threads '-p ${THREADS}'. It must be a positive integer."
if [[ -n "${RUN[download-wf]:-}" ]]; then
    if [[ ! "${DOWNLOAD_THREADS}" =~ ^[0-9]+$ ]] || (( DOWNLOAD_THREADS < 1 || DOWNLOAD_THREADS > 30 )); then
        cli_error "Invalid -j/--download-threads '${DOWNLOAD_THREADS}'. It must be an integer from 1 to 30." download-wf
    fi
    [[ "${MIN_CONTIG_LENGTH}" =~ ^[0-9]+$ ]] || cli_error "Invalid --min-contig-length '${MIN_CONTIG_LENGTH}'. It must be an integer." download-wf
    for (( I = 0; I < ${#QC_ARGS[@]}; I += 2 )); do
        case "${QC_ARGS[$I]}" in
            --max-contigs|--max-n-bases|--min-n50)
                [[ "${QC_ARGS[$((I + 1))]}" =~ ^[0-9]+$ ]] || cli_error "Invalid ${QC_ARGS[$I]} '${QC_ARGS[$((I + 1))]}'. It must be an integer." download-wf ;;
            *)
                if [[ ! "${QC_ARGS[$((I + 1))]}" =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v v="${QC_ARGS[$((I + 1))]}" 'BEGIN { exit !(v + 0 <= 100) }'; then
                    cli_error "Invalid ${QC_ARGS[$I]} '${QC_ARGS[$((I + 1))]}'. It must be a percentage from 0 to 100." download-wf
                fi ;;
        esac
    done
    [[ "${DEREP_MODE}" == "greedy" || "${DEREP_MODE}" == "dynamic" ]] \
        || cli_error "Invalid --derep-mode '${DEREP_MODE}'. Must be greedy or dynamic." download-wf
    if [[ ! "${DEREP_ANI}" =~ ^[0-9]+([.][0-9]+)?$ ]] || ! awk -v v="${DEREP_ANI}" 'BEGIN { exit !(v + 0 > 0 && v + 0 <= 100) }'; then
        cli_error "Invalid --derep-ani '${DEREP_ANI}'. It must be a percentage greater than 0 and at most 100 (e.g., 99.5)." download-wf
    fi
fi

# The selection query is checked now, before anything is run (the search module has the rules of every filter)
make_temp_dir WORK_DIR
QUERY_TMP="${WORK_DIR}/genome-selection-query.txt"
if [[ -n "${RUN[search-genomes]:-}" ]]; then
    {
        echo "# Prokfinder genome selection query (${PROKFINDER_VERSION})"
        echo "# Command: ${COMMAND_LINE}"
        echo "# FILTER <TAB> <id> <TAB> <key> <TAB> <value>: lines with the same id are one filter (its values combined with OR)"
        printf '%s\n' "${QUERY_LINES[@]}"
        printf 'EXPRESSION\t%s\n' "${EXPR_TOKENS[*]}"
    } > "${QUERY_TMP}"
    # (the taxon names of -N are checked with the NCBI taxonomy of the reference folder, when it is installed)
    bash "${MODULES_DIR}/search-genomes-in-reference-db.sh" -q "${QUERY_TMP}" --validate-only \
        --taxdump-dir "$(pf_ncbi_taxdump_dir "${NCBI_DIR}")" \
        || { echo "Use ${PROG} search-genomes -h (or --help) to see the options." >&2; exit 1; }
fi

# search-genomes and download-wf need the reference databases installed by setup
if [[ -z "${RUN[setup]:-}" && ( -n "${RUN[search-genomes]:-}" || -n "${RUN[download-wf]:-}" ) ]]; then
    if [[ "${#REFDB_MISSING[@]}" -gt 0 ]]; then
        cli_error "Reference databases not found in ${REF_DIR} (missing: $(for DB in NCBI BACDIVE GTDB; do [[ -n "${REFDB_MISSING[${DB}]:-}" ]] && printf "%s " "${DB}"; done | sed "s/ $//")). Install them first with: ${PROG} setup (it can be given together with the other sub-commands, e.g., ${PROG} setup search-genomes ...)."
    fi
    if [[ "${#REFDB_OUTDATED[@]}" -gt 0 ]]; then
        cli_error "The reference databases in ${REF_DIR} were installed by an older version of Prokfinder (tables to update: $(for DB in NCBI BACDIVE GTDB; do [[ -n "${REFDB_OUTDATED[${DB}]:-}" ]] && printf "%s " "${DB}"; done | sed "s/ $//")). Update them with: ${PROG} setup (the tables are rebuilt from the files already downloaded; nothing is downloaded again)."
    fi
fi

# download-wf works on the genomes selected by search-genomes (in this command or before, in the same results folder)
if [[ -n "${RUN[download-wf]:-}" && -z "${RUN[search-genomes]:-}" && ! -f "${SEARCH_CHECKPOINT}" ]]; then
    cli_error "download-wf needs the genomes selected by search-genomes, but ${SEARCH_CHECKPOINT} was not found. Run search-genomes first (it can be in the same command: ${PROG} search-genomes ... download-wf ...), with the same -o/--results-dir." download-wf
fi

############ EXECUTION ON THE CLUSTER (PBS) ############
#   With run-on-cluster the run is not made here: the command line without run-on-cluster is written into a PBS job
#   script and submitted with qsub. The job then runs this script on a compute node.

if [[ -n "${RUN[run-on-cluster]:-}" ]]; then
    PBS_STATUS=0
    if [[ "${PBS_ARGS[*]:-}" == *--test* ]]; then
        bash "${MODULES_DIR}/submit-prokfinder-pbs.sh" ${PBS_ARGS[@]+"${PBS_ARGS[@]}"} || PBS_STATUS=$?
    else
        bash "${MODULES_DIR}/submit-prokfinder-pbs.sh" ${PBS_ARGS[@]+"${PBS_ARGS[@]}"} -- "${JOB_ARGS[@]}" || PBS_STATUS=$?
    fi
    exit "${PBS_STATUS}"
fi

############ DIRECTORY STRUCTURE ############

OUTPUT_DIR="${PROKFINDER_DIR}/OUTPUT"
QUERY_FILE="${OUTPUT_DIR}/genome-selection-query.txt"
METADATA_DIR="${OUTPUT_DIR}/GENOME-METADATA"
METADATA="${METADATA_DIR}/genome-metadata.tsv"
GENOME_FILES_DIR="${OUTPUT_DIR}/GENOME-FILES"
RAW_DIR="${GENOME_FILES_DIR}/raw_genomes"
CLEAN_DIR="${GENOME_FILES_DIR}/clean_genomes"
QUALITY_DIR="${OUTPUT_DIR}/GENOME-QUALITY"
TAXONOMY_DIR="${OUTPUT_DIR}/GENOME-TAXONOMY"
HIERARCHY_DIR="${OUTPUT_DIR}/GENOME-DIR-HIERARCHY"
FINAL_DIR="${OUTPUT_DIR}/FINAL-RESULTS"

mkdir -p "${PROKFINDER_DIR}"
LOG_FILE="$(absolute_path "${PROKFINDER_DIR}")/prokfinder.log"
start_logging "${LOG_FILE}"
export PROKFINDER_LOG="${LOG_FILE}"

log "message" "================ Prokfinder ${PROKFINDER_VERSION} ================"
log "message" "Command: ${COMMAND_LINE}"
log "message" "Sub-commands: $(for S in setup search-genomes download-wf; do [[ -n "${RUN[${S}]:-}" ]] && printf '%s ' "${S}"; done)"
log "message" "Results folder: $(absolute_path "${PROKFINDER_DIR}") (log: ${LOG_FILE})"
mkdir -p "${REF_DIR}"
log "message" "Reference databases folder: $(absolute_path "${REF_DIR}")"

############ CHECKING THE DEPENDENCIES ############

# Reference databases whose module must run (setup): not installed, installed by an older Prokfinder, or forced; and
# NCBI when --ncbi-complete-metadata asks for a complete table that is missing
declare -A REFDB_NEEDED=()
refdb_needed() {
    local DB DOMAIN
    REFDB_NEEDED=()
    check_refdb_state
    for DB in NCBI BACDIVE GTDB; do
        if [[ "${FORCE_REFDB}" -eq 1 || -n "${REFDB_MISSING[${DB}]:-}" || -n "${REFDB_OUTDATED[${DB}]:-}" ]]; then
            REFDB_NEEDED["${DB}"]=1
        fi
    done
    if [[ "${NCBI_COMPLETE}" -eq 1 ]]; then
        for DOMAIN in "${PF_DOMAINS[@]}"; do
            [[ -s "$(pf_ncbi_complete_table "${NCBI_DIR}" "${DOMAIN}")" ]] || REFDB_NEEDED[NCBI]=1
        done
    fi
}
if [[ -n "${RUN[setup]:-}" ]]; then
    refdb_needed
fi

# Only the programs of the sub-commands given are required (the genome search does not need the programs of
# download-wf)
DEPENDENCY_ERRORS=0
check_program() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log "error" "'$1' was not found in the PATH (install: $(install_hint "$1")${2:+; $2})."
        DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
    fi
}

for TOOL in awk sort gzip curl; do
    check_program "${TOOL}"
done
if [[ -n "${REFDB_NEEDED[NCBI]:-}" ]]; then
    for TOOL in datasets dataformat taxonkit unzip tar; do
        check_program "${TOOL}"
    done
fi
if [[ -n "${REFDB_NEEDED[BACDIVE]:-}" ]]; then
    for TOOL in jq gawk taxonkit; do
        check_program "${TOOL}"
    done
fi
if [[ -n "${REFDB_NEEDED[GTDB]:-}" ]]; then
    check_program md5sum
fi
if [[ -n "${RUN[search-genomes]:-}" ]] && printf '%s\n' "${QUERY_LINES[@]}" | grep -q $'\tTAXON_NAME\t'; then
    check_program taxonkit "needed by -N/--ncbi-taxon-file"
fi
if [[ -n "${RUN[download-wf]:-}" ]]; then
    for TOOL in datasets seqkit checkm2 skder; do
        check_program "${TOOL}"
    done
    PROBLEM="$(checkm2db_problem)"
    if [[ -n "${PROBLEM}" ]]; then
        log "error" "${PROBLEM}"
        DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
    fi
    if [[ "${USE_GTDBTK}" -eq 1 ]]; then
        check_program gtdbtk
        PROBLEM="$(gtdbtk_data_problem)"
        if [[ -n "${PROBLEM}" ]]; then
            log "error" "${PROBLEM}"
            DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
        fi
    fi
fi

if [[ "${DEPENDENCY_ERRORS}" -gt 0 ]]; then
    die "${DEPENDENCY_ERRORS} problem(s) found in the dependencies. Fix them and run Prokfinder again (${PROG} --check lists every dependency)."
fi

log "message" "Dependencies:"
DEPS_LOGGED=("${DEPS_SEARCH[@]}")
if [[ -n "${RUN[download-wf]:-}" ]]; then
    DEPS_LOGGED+=("${DEPS_DOWNLOAD[@]}" "${DEPS_OPTIONAL[@]}")
fi
for ENTRY in "${DEPS_LOGGED[@]}"; do
    TOOL="${ENTRY%%|*}"
    if command -v "${TOOL}" >/dev/null 2>&1; then
        log "message" "  ${TOOL}: $(tool_version "${TOOL}")"
    fi
done
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    log "message" "  NCBI_API_KEY: set ('datasets' will use it: up to 10 requests per second)."
else
    log "message" "  NCBI_API_KEY: not set (optional; 'datasets' is limited to 5 requests per second)."
fi

# Runs a module; stops Prokfinder if it fails
run_module() {
    local NAME="$1"
    local RC=0
    shift
    log "message" "-------- ${NAME} --------"
    bash "${MODULES_DIR}/${NAME}" "$@" || RC=$?
    if [[ "${RC}" -ne 0 ]]; then
        die "Module ${NAME} failed (exit code ${RC}). Fix the problem shown above and run Prokfinder again: finished work is reused."
    fi
}

############ SETUP: REFERENCE DATABASES ############
#   The reference folder is shared by every run, so only one run installs it at a time: a run that finds another one
#   installing waits for it and then reuses what it installed.
#   With --force-refdb each database is installed in a temporary folder (<refdb>/.<DB>.new) and replaces the current
#   one only when its module finished without error: a failed download never leaves the folder without a database.
#   The databases already replaced are listed in <refdb>/.force-refdb.state, so the same command run again after a
#   failure skips them and goes on from the database that failed; the file is removed when the three are done.

FORCE_STATE="${REF_DIR}/.force-refdb.state"

# Runs the module of one reference database (into its folder, or into the temporary folder with --force-refdb)
# Usage: install_refdb <DB> <db_dir> <module> [module options]
install_refdb() {
    local DB="$1" DB_DIR="$2" MODULE="$3"
    local NEW_DIR="${REF_DIR}/.${DB}.new" OLD_DIR="${REF_DIR}/.${DB}.old"
    shift 3
    if [[ "${FORCE_REFDB}" -eq 0 ]]; then
        run_module "${MODULE}" -o "${DB_DIR}" "$@"
        return 0
    fi

    if grep -qx "${DB}" "${FORCE_STATE}" 2>/dev/null; then
        log "message" "${DB} was already installed again by the interrupted --force-refdb run (${FORCE_STATE}); skipping it."
        return 0
    fi

    if [[ -d "${NEW_DIR}" ]]; then
        log "message" "Resuming the new ${DB} installation in ${NEW_DIR} (left by an interrupted --force-refdb run)."
    else
        log "message" "Installing the new ${DB} database in the temporary folder ${NEW_DIR}; the current one is kept until it finishes."
    fi
    mkdir -p "${NEW_DIR}"
    run_module "${MODULE}" -o "${NEW_DIR}" "$@"

    # The new database is complete: it replaces the current one
    rm -rf "${OLD_DIR}"
    if [[ -d "${DB_DIR}" ]]; then
        mv "${DB_DIR}" "${OLD_DIR}"
    fi
    mv "${NEW_DIR}" "${DB_DIR}"
    sed -i "s#${NEW_DIR}#${DB_DIR}#g" "$(pf_completed_marker "${DB_DIR}")"
    rm -rf "${OLD_DIR}"
    echo "${DB}" >> "${FORCE_STATE}"
    log "message" "New ${DB} reference database installed in ${DB_DIR} (the previous one was removed)."
}

if [[ -n "${RUN[setup]:-}" ]]; then
    log "message" "======== setup: reference databases ========"
    if [[ "${FORCE_REFDB}" -eq 1 ]]; then
        log "message" "Option --force-refdb: the reference databases will be downloaded and installed again, each one in a temporary folder that replaces the current database only at the end."
    fi
    for DB in NCBI BACDIVE GTDB; do
        mkdir -p "${REF_DIR}/${DB}"
    done

    REFDB_LOCKED=0
    if [[ "${#REFDB_NEEDED[@]}" -gt 0 ]]; then
        if command -v flock >/dev/null 2>&1; then
            exec 9> "${REF_DIR}/.download.lock"
            if ! flock -n 9; then
                log "message" "Another Prokfinder run is installing the reference databases in ${REF_DIR}; waiting for it to finish..."
                flock 9
                [[ "${FORCE_REFDB}" -eq 1 ]] || refdb_needed
            fi
            REFDB_LOCKED=1
        else
            log "warning" "'flock' not found: two Prokfinder runs installing the reference databases at the same time could conflict."
        fi
    fi

    for DB in NCBI BACDIVE GTDB; do
        if [[ "${FORCE_REFDB}" -eq 0 && -n "${REFDB_OUTDATED[${DB}]:-}" ]]; then
            log "message" "The ${DB} reference tables were installed by an older version of Prokfinder: they are updated now from the files already downloaded (nothing is downloaded again)."
        fi
    done

    if [[ -n "${REFDB_NEEDED[NCBI]:-}" ]]; then
        NCBI_ARGS=(--assembly-source all)
        [[ "${NCBI_COMPLETE}" -eq 0 ]] || NCBI_ARGS+=(--complete-table)
        install_refdb NCBI "${NCBI_DIR}" download-reference-db-ncbi.sh "${NCBI_ARGS[@]}"
    else
        log "message" "NCBI reference data already installed ($(head -n 1 "$(pf_completed_marker "${NCBI_DIR}")"))."
    fi

    # The rank of the BacDive TaxIDs is found with the NCBI taxonomy downloaded with the NCBI tables. Tables of an
    # older Prokfinder are rebuilt from the strains already downloaded (no request to the BacDive API)
    if [[ -n "${REFDB_NEEDED[BACDIVE]:-}" ]]; then
        BACDIVE_ARGS=(--taxdump-dir "$(pf_ncbi_taxdump_dir "${NCBI_DIR}")")
        [[ "${FORCE_REFDB}" -eq 1 || -z "${REFDB_OUTDATED[BACDIVE]:-}" ]] || BACDIVE_ARGS+=(--skip-download)
        install_refdb BACDIVE "${BACDIVE_DIR}" download-reference-db-bacdive.sh "${BACDIVE_ARGS[@]}"
    else
        log "message" "BacDive reference data already installed ($(head -n 1 "$(pf_completed_marker "${BACDIVE_DIR}")"))."
    fi

    if [[ -n "${REFDB_NEEDED[GTDB]:-}" ]]; then
        install_refdb GTDB "${GTDB_DIR}" download-reference-db-gtdb.sh
    else
        log "message" "GTDB reference data already installed ($(head -n 1 "$(pf_completed_marker "${GTDB_DIR}")"))."
    fi
    # Every database was installed again: the next --force-refdb starts from scratch
    rm -f "${FORCE_STATE}"

    if [[ "${REFDB_LOCKED}" -eq 1 ]]; then
        flock -u 9
        exec 9>&-
    fi
    log "message" "Reference databases installed in $(absolute_path "${REF_DIR}")."
fi

############ SEARCH-GENOMES ############
#   The selection is reused when the query and the reference databases are the same as in the last search of this
#   results folder (recorded in SEARCH-GENOMES-COMPLETED.txt); --force-selection searches again.

if [[ -n "${RUN[search-genomes]:-}" ]]; then
    log "message" "======== search-genomes: genome selection ========"
    mkdir -p "${METADATA_DIR}"
    FINGERPRINT="$(
        grep -v '^#' "${QUERY_TMP}"
        for DB_DIR in "${NCBI_DIR}" "${BACDIVE_DIR}" "${GTDB_DIR}"; do
            echo "${DB_DIR}: $(file_mtime "$(pf_completed_marker "${DB_DIR}")")"
        done
    )"
    FINGERPRINT="$(printf '%s\n' "${FINGERPRINT}" | md5sum | cut -d' ' -f1)"

    if [[ "${FORCE_SELECTION}" -eq 1 ]]; then
        log "message" "Option --force-selection: removing ${SEARCH_CHECKPOINT} and searching again."
        rm -f "${SEARCH_CHECKPOINT}"
    fi

    if [[ -s "${METADATA}" && -f "${SEARCH_CHECKPOINT}" ]] && grep -qx "fingerprint=${FINGERPRINT}" "${SEARCH_CHECKPOINT}"; then
        log "message" "The selection filters and reference databases did not change since the last search: reusing ${METADATA} (use --force-selection to search again)."
    else
        rm -f "${SEARCH_CHECKPOINT}"
        log "message" "Selection query:"
        while IFS= read -r LINE; do
            log "message" "    ${LINE//$'\t'/ }"
        done < <(grep -v '^#' "${QUERY_TMP}")

        run_module search-genomes-in-reference-db.sh -q "${QUERY_TMP}" -r "${REF_DIR}" -o "${METADATA}"
        cp "${QUERY_TMP}" "${QUERY_FILE}"
        {
            echo "search-genomes completed on $(date +"%Y-%m-%d %H:%M:%S")"
            echo "command=${COMMAND_LINE}"
            echo "expression=${EXPR_TOKENS[*]}"
            echo "assembly_source=${ASSEMBLY_SOURCE}"
            echo "selected_genomes=$(count_records "${METADATA}")"
            echo "metadata=${METADATA}"
            echo "fingerprint=${FINGERPRINT}"
        } > "${SEARCH_CHECKPOINT}"
    fi
    log "message" "Selected genome assemblies: $(count_records "${METADATA}") (metadata: ${METADATA})"
fi

if [[ -z "${RUN[download-wf]:-}" ]]; then
    log "message" "Prokfinder finished successfully."
    exit 0
fi

############ DOWNLOAD-WF ############

log "message" "======== download-wf: download, quality control, taxonomy and dereplication ========"
[[ -s "${METADATA}" ]] || die "The table of the selected genomes was not found: ${METADATA}. Run search-genomes again."
ASSEMBLY_SOURCE="$(sed -n 's/^assembly_source=//p' "${SEARCH_CHECKPOINT}" | head -n 1)"
ASSEMBLY_SOURCE="${ASSEMBLY_SOURCE:-genbank}"

SELECTED_LIST="${METADATA_DIR}/selected-accessions.txt"
tail -n +2 "${METADATA}" | cut -f1 > "${SELECTED_LIST}"
N_SELECTED="$(count_lines "${SELECTED_LIST}")"
log "message" "Genomes selected by search-genomes: ${N_SELECTED} ($(sed -n 's/^search-genomes //p' "${SEARCH_CHECKPOINT}" | head -n 1))"

if [[ "${N_SELECTED}" -eq 0 ]]; then
    log "message" "No genome assembly was selected by search-genomes, so there is nothing to download."
    exit 0
fi

if [[ "${FORCE_CURATION}" -eq 1 ]]; then
    log "message" "Option --force-curation: removing the clean genomes and the quality, taxonomy, hierarchy and final results."
    rm -rf "${CLEAN_DIR}" "${QUALITY_DIR}" "${TAXONOMY_DIR}" "${HIERARCHY_DIR}" "${FINAL_DIR}"
fi
mkdir -p "${RAW_DIR}" "${CLEAN_DIR}" "${QUALITY_DIR}" "${TAXONOMY_DIR}" "${HIERARCHY_DIR}" \
    "${FINAL_DIR}/STATISTICS/TABLES" "${FINAL_DIR}/STATISTICS/GRAPHS"

GZIP_FLAG=()
if [[ "${GZIP}" -eq 1 ]]; then
    GZIP_FLAG=(-z)
fi

# Thresholds drawn as reference lines in the graphs (plot-statistics.R): the values of this run, never inherited
# from the environment
export QC_MIN_COMPLETENESS=75 QC_MAX_CONTAMINATION=5 QC_MIN_N50=5000 QC_MAX_CONTIGS=1000
for (( I = 0; I < ${#QC_ARGS[@]}; I += 2 )); do
    case "${QC_ARGS[$I]}" in
        --min-completeness)  QC_MIN_COMPLETENESS="${QC_ARGS[$((I + 1))]}" ;;
        --max-contamination) QC_MAX_CONTAMINATION="${QC_ARGS[$((I + 1))]}" ;;
        --min-n50)           QC_MIN_N50="${QC_ARGS[$((I + 1))]}" ;;
        --max-contigs)       QC_MAX_CONTIGS="${QC_ARGS[$((I + 1))]}" ;;
    esac
done

# Writes the accessions of the metadata rows that satisfy an awk condition on the named columns.
# Usage: accessions_where <output_file> '<condition with $COL_NAME>' <COL_NAME> [COL_NAME ...]
#        e.g. accessions_where list.txt '$QC_STATUS == "PASSED"' QC_STATUS
accessions_where() {
    local OUT="$1"
    local CONDITION="$2"
    local NAME INDEX
    shift 2
    for NAME in "$@"; do
        INDEX="$(require_column "${METADATA}" "${NAME}")"
        CONDITION="${CONDITION//\$${NAME}/\$${INDEX}}"
    done
    awk -F'\t' "NR > 1 && (${CONDITION}) { print \$1 }" "${METADATA}" > "${OUT}"
}

# Runs the final summary and stops (used when a step leaves no genome for the next ones)
finish_early() {
    log "warning" "$1"
    run_module summarize-final-results.sh -m "${METADATA}" -r "${NCBI_DIR}" -s "${ASSEMBLY_SOURCE}" -o "${FINAL_DIR}"
    exit 0
}

# Step 2: genome download
run_module download-selected-genomes-from-db.sh -l "${SELECTED_LIST}" -r "${NCBI_DIR}" -o "${RAW_DIR}" \
    -t "${GENOME_FILES_DIR}/genome-download-status.tsv" -j "${DOWNLOAD_THREADS}" ${GZIP_FLAG[@]+"${GZIP_FLAG[@]}"}
add_metadata_columns "${METADATA}" "${GENOME_FILES_DIR}/genome-download-status.tsv"

DOWNLOADED_LIST="${WORK_DIR}/downloaded.txt"
accessions_where "${DOWNLOADED_LIST}" '$DOWNLOAD_STATUS == "DOWNLOADED"' DOWNLOAD_STATUS
if [[ ! -s "${DOWNLOADED_LIST}" ]]; then
    finish_early "No selected genome could be downloaded."
fi

# Step 3: contig filtering
run_module filter-low-quality-contigs-from-genomes.sh -l "${DOWNLOADED_LIST}" -i "${RAW_DIR}" -o "${CLEAN_DIR}" \
    -t "${QUALITY_DIR}/contig-filtering.tsv" -m "${MIN_CONTIG_LENGTH}" -p "${THREADS}" ${GZIP_FLAG[@]+"${GZIP_FLAG[@]}"}
add_metadata_columns "${METADATA}" "${QUALITY_DIR}/contig-filtering.tsv"

CLEAN_LIST="${WORK_DIR}/clean.txt"
accessions_where "${CLEAN_LIST}" '$CONTIG_FILTER_STATUS == "OK"' CONTIG_FILTER_STATUS
if [[ ! -s "${CLEAN_LIST}" ]]; then
    finish_early "No genome has contigs left after the contig filtering."
fi

# Step 4: quality control
run_module evaluate-and-filter-genome-by-qc.sh -l "${CLEAN_LIST}" -i "${CLEAN_DIR}" -o "${QUALITY_DIR}" -p "${THREADS}" \
    ${QC_ARGS[@]+"${QC_ARGS[@]}"}
add_metadata_columns "${METADATA}" "${QUALITY_DIR}/genome-quality.tsv"

HIGH_QUALITY_LIST="${WORK_DIR}/high_quality.txt"
accessions_where "${HIGH_QUALITY_LIST}" '$QC_STATUS == "PASSED"' QC_STATUS
if [[ ! -s "${HIGH_QUALITY_LIST}" ]]; then
    finish_early "No genome passed the quality filters."
fi

# Step 5: GTDB taxonomy
GTDBTK_FLAG=()
if [[ "${USE_GTDBTK}" -eq 1 ]]; then
    GTDBTK_FLAG=(--gtdbtk)
fi
run_module assign-taxonomy-gtdb.sh -l "${HIGH_QUALITY_LIST}" -g "${GTDB_DIR}" -o "${TAXONOMY_DIR}" -i "${CLEAN_DIR}" \
    -p "${THREADS}" ${GTDBTK_FLAG[@]+"${GTDBTK_FLAG[@]}"}
add_metadata_columns "${METADATA}" "${TAXONOMY_DIR}/genome-taxonomy.tsv"

CLASSIFIED_LIST="${WORK_DIR}/classified.txt"
accessions_where "${CLASSIFIED_LIST}" '$QC_STATUS == "PASSED" && $GTDB_TAXONOMY_SOURCE != "NULL"' QC_STATUS GTDB_TAXONOMY_SOURCE
if [[ ! -s "${CLASSIFIED_LIST}" ]]; then
    finish_early "No high-quality genome has a GTDB taxonomy."
fi

# Step 6: taxonomic directory hierarchy
run_module build-taxonomic-dir-hierarchy.sh -m "${METADATA}" -l "${CLASSIFIED_LIST}" -o "${HIERARCHY_DIR}"
add_metadata_columns "${METADATA}" "${HIERARCHY_DIR}/genome-dir-hierarchy.tsv"

# Step 7: dereplication
MGE_FLAG=()
if [[ "${FILTER_MGE}" -eq 1 ]]; then
    MGE_FLAG=(--filter-mge)
fi
run_module dereplicate-genomes-skder.sh -d "${HIERARCHY_DIR}" -m "${METADATA}" \
    -t "${HIERARCHY_DIR}/genome-dereplication.tsv" -p "${THREADS}" --mode "${DEREP_MODE}" --ani "${DEREP_ANI}" ${MGE_FLAG[@]+"${MGE_FLAG[@]}"}
add_metadata_columns "${METADATA}" "${HIERARCHY_DIR}/genome-dereplication.tsv"

# Step 8: final results
run_module summarize-final-results.sh -m "${METADATA}" -r "${NCBI_DIR}" -s "${ASSEMBLY_SOURCE}" -o "${FINAL_DIR}"

log "message" "Prokfinder finished successfully. Genome metadata (every step): ${METADATA}"
log "message" "High-quality representative genomes: ${FINAL_DIR}/high-quality-selected-genomes.txt"
