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

PROKFINDER_VERSION="0.3.0"
PROG="$(basename "$0")"

############ HELP ############

help_global_options() {
    cat << EOF
Global options (accepted anywhere in the command line, with any sub-command):
  -o, --results-dir <dir>   Results folder of this run (default: ./prokfinder-results)
  --refdb-dir <dir>         Folder of the reference databases, shared by every run (default: \$HOME/prokfinder-refdb)
  -p, --threads <n>         CPUs used by every program: seqkit, CheckM2, GTDB-Tk, skDER (default: 4; with
                            run-on-cluster, the CPUs of the job)
  -h, --help                Show the help (of the sub-command, when given after it)
EOF
}

help_main() {
    cat << EOF
Prokfinder ${PROKFINDER_VERSION} - selection, download and curation of prokaryote genomes

Usage: ${PROG} [GLOBAL OPTIONS] <sub-command> [options] [<sub-command> [options] ...]

Sub-commands (when several are given, they run in this order):
  setup            Install the reference databases (NCBI, BacDive and GTDB). Required once, before search-genomes
                   and download-wf; also updates databases installed by an older Prokfinder
  search-genomes   Select the genomes that match the selection filters (NCBI, BacDive and GTDB, joined with AND/OR)
  download-wf      Download the selected genomes, remove small contigs, filter them by quality (seqkit, CheckM2),
                   assign the GTDB taxonomy and dereplicate each taxon (skDER). Needs a finished search-genomes
  run-on-cluster   Run the other sub-commands as a job on a compute node of a cluster (PBS)

$(help_global_options)

Examples:
  ${PROG} setup
  ${PROG} search-genomes -T 45:80 AND -O rod,coccus
  ${PROG} -o thermophiles setup search-genomes -T 45:80 AND -O rod,coccus download-wf --min-completeness 90
  ${PROG} run-on-cluster --queue qtime --walltime 24:00:00 search-genomes -N taxa.txt OR -i "hot spring"

Help of a sub-command: ${PROG} <sub-command> -h
EOF
}

help_setup() {
    cat << EOF
Usage: ${PROG} [GLOBAL OPTIONS] setup [--force-refdb] [<other sub-commands>]

Installs the reference databases in the folder of --refdb-dir (default: \$HOME/prokfinder-refdb):
  NCBI     metadata of every Bacteria/Archaea assembly (GenBank and RefSeq) and the NCBI taxonomy
  BacDive  phenotypic traits of the strains, merged by strain and by NCBI TaxID
  GTDB     taxonomy and metadata of the genomes classified by GTDB
The download takes hours and is made only once: the databases are shared by every run and every results folder.
Databases already installed are not downloaded again; databases installed by an older Prokfinder are updated from
the files already downloaded (nothing is downloaded again).

Options:
  --force-refdb   Download and install the reference databases again. Each one is installed in a temporary folder
                  (<refdb>/.<DB>.new) and replaces the current database only when it finished without error; if
                  it fails, the current database is kept and the same command resumes the new installation
                  (the databases already replaced are not downloaded again)

$(help_global_options)
EOF
}

help_search() {
    cat << EOF
Usage: ${PROG} [GLOBAL OPTIONS] search-genomes <filter> [AND|OR <filter> ...] [options]

Selects, in the NCBI tables, the genomes that match the selection filters. Writes
OUTPUT/GENOME-METADATA/genome-metadata.tsv and the checkpoint SEARCH-GENOMES-COMPLETED.txt in the results folder.

How the filters are written:
  1. Each flag below is one filter. At least one filter is required.
  2. Several values of the same flag are separated by commas, without spaces, and are combined with OR: the filter
     matches when ANY of the values matches (e.g., -t 562,1280). Each value that contains spaces is written between
     its own double quotes (e.g., -i "hot spring","hydrothermal vent",geyser).
  3. Between two filters, AND or OR is required:
       -i "hot spring",geyser AND -c Brazil   genomes from a hot spring or a geyser, and from Brazil
       -i "hot spring",geyser OR -c Brazil    genomes from a hot spring or a geyser, or from Brazil
  4. AND is evaluated before OR. Parentheses group filters; they must be quoted or escaped for the shell:
       -T 45:80 AND "(" -O rod OR -M anaerobe ")"
  5. List files (-f, -A, -N, -I, -C): one value per line; spaces at the start/end of the lines and Windows line
     endings (\r) are removed, and empty lines are skipped.

NCBI filters:
  -t, --taxid <id>[,<id>]          NCBI TaxID of the organism or of any rank of its lineage
                                   (e.g., 1239 selects every Bacillota genome)
  -f, --taxid-file <file>          NCBI TaxIDs, one per line
  -a, --accession <acc>[,<acc>]    Accession with version, GenBank (GCA_) or RefSeq (GCF_) (e.g., GCA_000005845.2)
  -A, --accession-file <file>      Accessions, one per line; GCA_ and GCF_ may be mixed (GCF_x.1 selects its GenBank
                                   pair GCA_x.1). Every line that is not an accession is reported
  --ncbi-taxon <name>[,<name>]     Taxon name of the NCBI lineage, any rank (e.g., Thermus)
  -N, --ncbi-taxon-file <file>     NCBI taxon names, one per line (e.g., Thermus aquaticus), converted to TaxIDs with
                                   taxonkit name2taxid; a name that is not a Bacteria/Archaea taxon is an error
  -i, --isolation-source <text>[,<text>]
                                   Text contained in the isolation source (e.g., "hot spring","hydrothermal vent")
  -I, --isolation-source-file <f>  Isolation sources, one per line
  -c, --country <name>[,<name>]    Country/place of the geographic location (e.g., Brazil,Iceland)
  -C, --country-file <file>        Countries, one per line
  -g, --genome-type <type>         isolate_genome or mag
  -G, --gc <min:max>               GC content (%) of the assembly (e.g., 40:65; either side may be empty)

BacDive filters (phenotype; ranges <min:max>, either side may be empty):
  -T, --temp-opt <min:max>         Optimal growth temperature in °C (e.g., 45:80)
  -P, --ph-opt <min:max>           Optimal pH (e.g., 5.5:8.5)
  -S, --salt-opt <min:max>         Optimal NaCl concentration in % (e.g., 0:5)
  --temp-growth <min:max>          Temperature at which growth was recorded (most strains only have this one)
  --ph-growth <min:max>            pH at which growth was recorded
  --salt-growth <min:max>          NaCl concentration at which growth was recorded
  -O, --cell-shape <class>[,...]   Cell shape (e.g., rod,coccus)
  -M, --oxygen <class>[,...]       Oxygen tolerance (e.g., anaerobe)
  -B, --bacdive-class <trait=class[,class]>
                                   Any trait divided in classes (e.g., -B gram_stain=negative)
    Traits and classes: cell_shape (ovoid spiral coccus curved rod pleomorphic filament other),
    oxygen_tolerance (aerobe anaerobe facultative microaerophile), gram_stain (positive negative),
    motility / pigmentation_production / forms_multicellular_complex / spore_formation (yes no),
    flagellum_arrangement (polar peritrichous gliding), colony_color (colorless cream yellow white brown pink
    orange red grey green blue), energy_source (chemotroph phototroph), electron_donor (organotroph lithotroph),
    carbon_source (autotroph heterotroph), biosafety_level (1 2 3 4)
  BacDive is used with one row per NCBI TaxID: the references of each strain and then the strains of the same TaxID
  are merged (numbers: mean; classes: share of the references). A number matches when the mean is inside the range;
  a class matches when it was recorded for the TaxID. A TaxID of rank species selects the genomes of the species;
  a TaxID below the species (strain) selects the genomes of that strain; TaxIDs above the species are not used.

GTDB filters:
  --gtdb-taxon <name>[,<name>]     Taxon name of the GTDB lineage, any rank (e.g., g__Thermus or Thermus)

Options:
  -s, --assembly-source <src>      NCBI table searched: genbank (default) or refseq
  --force-selection                Search again (by default, the selection is reused when the filters and the
                                   reference databases did not change)

Result columns (genome-metadata.tsv), one per filter:
  SELECTED_BY_<FLAG>-<OPTION>  flags with a single option, TRUE/FALSE (e.g., SELECTED_BY_GENOME_TYPE-ISOLATE-GENOME)
  SELECTED_BY_<FLAG>-<RANGE>   numeric flags, TRUE/FALSE (e.g., SELECTED_BY_OPT_TEMP-45:80)
  SELECTED_BY_<FLAG>           flags with several values, TRUE/FALSE, plus a column with the values that matched,
                               separated by ";" (e.g., SELECTED_BY_ISOLATION_SOURCE and ISOLATION_SOURCE_KEYWORD_MATCH)

$(help_global_options)

Examples:
  ${PROG} search-genomes -T 45:80 AND -O rod,coccus AND -g isolate_genome
  ${PROG} search-genomes -N taxa.txt OR "(" --temp-growth 70: AND -i "hot spring","hydrothermal vent" ")"
  ${PROG} search-genomes --gtdb-taxon g__Thermus AND -c Brazil,Iceland -s refseq
EOF
}

help_download() {
    cat << EOF
Usage: ${PROG} [GLOBAL OPTIONS] download-wf [options]

Runs the workflow on the genomes selected by search-genomes (in this command or in a previous one with the same
results folder: it needs SEARCH-GENOMES-COMPLETED.txt in the results folder):
  download -> small contigs removed -> quality filters -> GTDB taxonomy -> directory hierarchy -> dereplication
  -> final tables and graphs

Download:
  -j, --download-threads <n>       Concurrent genome downloads, 1-30 (default: 10)
  -z, --gzip                       Keep the genome files compressed (.fna.gz)

Contig filtering:
  --min-contig-length <bp>         Contigs shorter than this are removed (default: 1000)

Genome quality filters (seqkit + CheckM2):
  --max-contigs <n>                The genome must have fewer contigs than this (default: 1000)
  --max-n-bases <n>                The genome must have fewer ambiguous bases (N) than this (default: 100000)
  --min-n50 <bp>                   Minimum N50 (default: 5000)
  --min-completeness <pct>         Minimum CheckM2 completeness (default: 75)
  --max-contamination <pct>        Maximum CheckM2 contamination (default: 5)

Taxonomy:
  --gtdbtk                         Classify with GTDB-Tk the genomes not found in the GTDB tables
                                   (WARNING: GTDB-Tk needs more than 200 GB of RAM); by default they are only listed
                                   in OUTPUT/GENOME-TAXONOMY/NON-CLASSIFIED-GENOMES-BY-GTDB.txt

Dereplication (skDER, inside each taxon; aligned fraction 90%):
  --derep-ani <pct>                ANI (% identity) above which two genomes are redundant (default: 99.5, the
                                   skDER default); lower values keep fewer representatives per taxon
  --derep-mode <mode>              greedy (default) or dynamic
  --filter-mge                     Mask the mobile genetic elements (prophages, plasmids; predicted with PhiSpy)
                                   before comparing the genomes

Re-running:
  --force-curation                 Discard the clean genomes and every quality/taxonomy/dereplication result (the
                                   downloaded genomes are kept). By default, finished work is reused

Programs: datasets, seqkit, checkm2 (CHECKM2DB must point to its database file), skder, gtdbtk (with --gtdbtk;
GTDBTK_DATA_PATH must point to its reference package). On the cluster, these are set by the environment modules.

$(help_global_options)

Examples:
  ${PROG} download-wf --min-completeness 90 --max-contamination 2
  ${PROG} search-genomes -T 45:80 AND -O rod download-wf --derep-ani 99 --derep-mode dynamic --filter-mge -z
EOF
}

help_cluster() {
    cat << EOF
Usage: ${PROG} [GLOBAL OPTIONS] run-on-cluster [options] <sub-commands>

Runs the other sub-commands of the command line as a job on a compute node, with the PBS queue manager: a .pbs
script that loads the environment modules of the cluster and runs Prokfinder is written to <results>/PBS-JOBS/ and
submitted with qsub. The command line is checked before the submission.

Options:
  --queue <queue>         Submission queue (default: roteamento):
                            qtime       tests and short runs (walltime of at most 72:00:00)
                            roteamento  common runs (at most 24024 MB of RAM per job)
                            high_mem    runs that need a lot of RAM (25 GB or more; e.g., download-wf --gtdbtk)
  --mem <ram>             RAM of the job: 20gb, 20000mb or only the number in GB (default: 20gb)
  --cpus <n>              CPUs of the job (default: 20); Prokfinder uses all of them unless -p is given
  --walltime <HH:MM:SS>   Maximum run time of the job (default: 2160:00:00; 72:00:00 with qtime)
  --test                  Submit only a short test job that loads the modules and shows the version of every
                          program (no other sub-command is needed)
  --dry-run               Only write the .pbs script (do not submit it)

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
GLOBAL_FLAGS=" -h --help "
GLOBAL_VALUES=" -o --results-dir --refdb-dir -p --threads "
SETUP_FLAGS=" --force-refdb "
SETUP_VALUES=" "
SEARCH_FLAGS=" --force-selection "
SEARCH_VALUES=" -s --assembly-source "
SEARCH_FILTERS=" -t --taxid -f --taxid-file -a --accession -A --accession-file --ncbi-taxon -N --ncbi-taxon-file
 -i --isolation-source -I --isolation-source-file -c --country -C --country-file -g --genome-type -G --gc
 -T --temp-opt -P --ph-opt -S --salt-opt --temp-growth --ph-growth --salt-growth -O --cell-shape -M --oxygen
 -B --bacdive-class --gtdb-taxon "
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
FORCE_REFDB=0
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
        # setup
        --force-refdb)                FORCE_REFDB=1 ;;
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
        -B|--bacdive-class)           start_filter "${ARG}"; add_value BACDIVE_CLASS "${VAL}" ;;
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

############ DEPENDENCIES ############

# Reference databases whose module must run (setup): not installed, installed by an older Prokfinder, or forced
declare -A REFDB_NEEDED=()
refdb_needed() {
    local DB
    REFDB_NEEDED=()
    check_refdb_state
    for DB in NCBI BACDIVE GTDB; do
        if [[ "${FORCE_REFDB}" -eq 1 || -n "${REFDB_MISSING[${DB}]:-}" || -n "${REFDB_OUTDATED[${DB}]:-}" ]]; then
            REFDB_NEEDED["${DB}"]=1
        fi
    done
}
if [[ -n "${RUN[setup]:-}" ]]; then
    refdb_needed
fi

# Prints the first line printed by a program to get its version
tool_version() {
    local VERSION
    VERSION="$("$@" 2>&1 | head -n 1)" || true
    echo "${VERSION:-(version not reported)}"
}

DEPENDENCY_ERRORS=0
check_program() {
    if ! command -v "$1" >/dev/null 2>&1; then
        log "error" "'$1' was not found in the PATH ($2)."
        DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
    fi
}

for TOOL in awk sort gzip curl; do
    check_program "${TOOL}" "basic Unix tool"
done
if [[ -n "${REFDB_NEEDED[NCBI]:-}" ]]; then
    check_program datasets "conda install -c conda-forge ncbi-datasets-cli"
    check_program dataformat "conda install -c conda-forge ncbi-datasets-cli"
    check_program taxonkit "conda install -c bioconda taxonkit"
    check_program unzip "apt install unzip"
fi
if [[ -n "${REFDB_NEEDED[BACDIVE]:-}" ]]; then
    check_program jq "conda install -c conda-forge jq"
    check_program gawk "conda install -c conda-forge gawk"
    check_program taxonkit "conda install -c bioconda taxonkit"
fi
if [[ -n "${REFDB_NEEDED[GTDB]:-}" ]]; then
    check_program md5sum "coreutils"
fi
if [[ -n "${RUN[search-genomes]:-}" ]] && printf '%s\n' "${QUERY_LINES[@]}" | grep -q $'\tTAXON_NAME\t'; then
    check_program taxonkit "conda install -c bioconda taxonkit (needed by -N/--ncbi-taxon-file)"
fi
if [[ -n "${RUN[download-wf]:-}" ]]; then
    check_program datasets "conda install -c conda-forge ncbi-datasets-cli"
    check_program seqkit "conda install -c bioconda seqkit"
    check_program checkm2 "conda install -c bioconda checkm2"
    check_program skder "conda install -c bioconda skder"

    if [[ -z "${CHECKM2DB:-}" ]]; then
        log "error" "CHECKM2DB is not set. It must be the path of the CheckM2 database file, e.g. export CHECKM2DB=/path/CheckM2_database/uniref100.KO.1.dmnd (download it with: checkm2 database --download; on the cluster, the checkm2 module sets it)."
        DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
    elif [[ ! -f "${CHECKM2DB}" || ! -r "${CHECKM2DB}" ]]; then
        log "error" "CHECKM2DB must point to the database file itself (uniref100.KO.1.dmnd); not found or not readable: ${CHECKM2DB}"
        DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
    fi

    if [[ "${USE_GTDBTK}" -eq 1 ]]; then
        check_program gtdbtk "conda install -c bioconda gtdbtk"
        if [[ -z "${GTDBTK_DATA_PATH:-}" || ! -d "${GTDBTK_DATA_PATH}" ]]; then
            log "error" "GTDBTK_DATA_PATH is not set or is not a folder. It must point to the GTDB-Tk reference package (on the cluster, the gtdbtk module sets it)."
            DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
        fi
    fi
fi

if [[ "${DEPENDENCY_ERRORS}" -gt 0 ]]; then
    die "${DEPENDENCY_ERRORS} problem(s) found in the dependencies. Fix them and run Prokfinder again."
fi

log "message" "Dependencies:"
for TOOL in datasets dataformat taxonkit seqkit checkm2 skder gtdbtk Rscript; do
    if command -v "${TOOL}" >/dev/null 2>&1; then
        case "${TOOL}" in
            seqkit)     log "message" "  seqkit: $(tool_version seqkit version)" ;;
            dataformat) log "message" "  dataformat: $(tool_version dataformat version)" ;;
            taxonkit)   log "message" "  taxonkit: $(tool_version taxonkit version)" ;;
            Rscript)    log "message" "  Rscript: $(tool_version Rscript --version)" ;;
            *)          log "message" "  ${TOOL}: $(tool_version "${TOOL}" --version)" ;;
        esac
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
        install_refdb NCBI "${NCBI_DIR}" download-reference-db-ncbi.sh --assembly-source all
    else
        log "message" "NCBI reference data already installed ($(head -n 1 "$(pf_completed_marker "${NCBI_DIR}")"))."
    fi

    # The rank of the BacDive TaxIDs is found with the NCBI taxonomy downloaded with the NCBI tables
    if [[ -n "${REFDB_NEEDED[BACDIVE]:-}" ]]; then
        install_refdb BACDIVE "${BACDIVE_DIR}" download-reference-db-bacdive.sh --taxdump-dir "$(pf_ncbi_taxdump_dir "${NCBI_DIR}")"
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
