#!/usr/bin/env bash
#
# main-prokfinder.sh - Prokfinder: pipeline for downloading, curating and identifying genomes of cultivated and
# uncultivated prokaryotes.
#
# This script only controls the execution: it reads the options, checks the programs and environment variables,
# creates the directory structure and calls the modules (folder modules/), each one doing one step:
#
#   0. download-reference-db-ncbi.sh / -bacdive.sh / -gtdb.sh   reference database (only on the first run)
#   1. search-genomes-in-reference-db.sh                        genome search (genome-selection-query.txt)
#   2. download-selected-genomes-from-db.sh                     download (rehydration) of the selected genomes
#   3. filter-low-quality-contigs-from-genomes.sh               contigs < 1000 bp removed, sorted and renamed
#   4. evaluate-and-filter-genome-by-qc.sh                      seqkit + CheckM2 quality filters
#   5. assign-taxonomy-gtdb.sh                                  GTDB taxonomy (GTDB-Tk optional)
#   6. build-taxonomic-dir-hierarchy.sh                         Domain/.../Species directory hierarchy
#   7. dereplicate-genomes-skder.sh                             dereplication inside each taxon (skDER)
#   8. summarize-final-results.sh                               final tables and graphs
#
# Every step adds its columns to a single table, OUTPUT/GENOME-METADATA/genome-metadata.tsv (one row per genome).
# Every step can be re-executed: finished work is reused (see the "Re-running" section of the help).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODULES_DIR="${SCRIPT_DIR}/modules"
# shellcheck source=lib/common.sh
source "${SCRIPT_DIR}/lib/common.sh"

PROKFINDER_VERSION="0.2.0"

usage() {
    cat << EOF
Prokfinder ${PROKFINDER_VERSION} - selection, download and curation of prokaryote genomes

Usage: $(basename "$0") [SELECTION FILTERS] [OPTIONS]

  Rules:
    1. At least one selection filter (NCBI, BacDive or GTDB) is required.
    2. To give several values to the same flag, separate them with commas, without spaces (e.g., -t 562,1280).
       A value that contains spaces must be written between double quotes (e.g., -i "hot spring,geyser").
    3. A genome is selected when it matches AT LEAST ONE of the values given to a flag (or listed in its file)
       and, at the same time, EVERY flag used.
       Example: -i "hot spring,geyser" -c Brazil selects genomes isolated from a hot spring or from a geyser,
       but only those from Brazil.

  NCBI filters:
    -t, --taxid <id>[,<id>]          NCBI TaxID of the organism or of any rank of its lineage
                                     (e.g., 1239 selects every Bacillota genome)
    -f, --taxid-file <file>          NCBI TaxIDs, one per line
    -a, --accession <acc>[,<acc>]    GenBank (GCA_) or RefSeq (GCF_) accession; the version is optional
    -A, --accession-file <file>      Accessions, one per line
    --ncbi-taxon <name>[,<name>]     Taxon name of the NCBI lineage, any rank (e.g., Thermus)
    -i, --isolation-source <text>[,<text>]
                                     Text contained in the isolation source (e.g., "hot spring,hydrothermal vent")
    -I, --isolation-source-file <f>  Isolation sources, one per line
    -c, --country <name>[,<name>]    Country/place of the geographic location (e.g., Brazil,Iceland)
    -C, --country-file <file>        Countries, one per line
    -g, --genome-type <type>         isolate_genome, mag or all (default: all)
    -G, --gc <min:max>               GC content (%) of the assembly (e.g., 40:65)
    --ncbi-range <COLUMN=min:max>    Range on any numeric column of the NCBI table
                                     (e.g., ASSEMBLY_STATS_CONTIG_N50=50000:)
    --ncbi-match <COLUMN=text>       Text contained in any column of the NCBI table
                                     (e.g., ASSEMBLY_LEVEL="Complete Genome")
    -s, --assembly-source <src>      NCBI table searched: genbank (default) or refseq

  BacDive filters (phenotype; linked to NCBI by the TaxID). Ranges: <min:max>, either side may be empty.
    -T, --temp-opt <min:max>         Optimal growth temperature in °C (e.g., 45:80)
    -P, --ph-opt <min:max>           Optimal pH (e.g., 5.5:8.5)
    -S, --salt-opt <min:max>         Optimal NaCl concentration in % (e.g., 0:5)
    --temp-growth <min:max>          Any temperature at which growth was recorded (minimum, maximum, optimum...)
    --ph-growth <min:max>            Any pH at which growth was recorded
    --salt-growth <min:max>          Any NaCl concentration at which growth was recorded
    --bacdive-range <column=min:max> Range on any numeric BacDive column (e.g., temperature_delta=20:)
    -O, --cell-shape <class>[,...]   Cell shape (e.g., rod,coccus)
    -M, --oxygen <class>[,...]       Oxygen tolerance (e.g., anaerobe)
    -B, --bacdive-class <trait=class[,class]>
                                     Any trait divided in classes; can be repeated for different traits
                                     (e.g., -B gram_stain=negative -B energy_source=chemotroph)
      Traits and classes: cell_shape (ovoid spiral coccus curved rod pleomorphic filament other),
      oxygen_tolerance (aerobe anaerobe facultative microaerophile), gram_stain (positive negative),
      motility / pigmentation_production / forms_multicellular_complex / spore_formation (yes no),
      flagellum_arrangement (polar peritrichous gliding), colony_color (colorless cream yellow white brown pink
      orange red grey green blue), energy_source (chemotroph phototroph), electron_donor (organotroph lithotroph),
      carbon_source (autotroph heterotroph), biosafety_level (1 2 3 4)
    A strain matches when at least ONE of its recorded values satisfies the filter (references may disagree);
    strains with no value for a filtered trait are excluded.

  GTDB filters:
    --gtdb-taxon <name>[,<name>]     Taxon name of the GTDB lineage, any rank (e.g., g__Thermus or Thermus)
    --gtdb-range <column=min:max>    Range on any numeric column of the GTDB metadata
                                     (e.g., checkm2_completeness=95:)
    --gtdb-match <column=text>       Text contained in any column of the GTDB metadata

  Execution:
    -o, --results-dir <dir>          Results folder of this run: selected genomes, quality, taxonomy and results
                                     (default: ./prokfinder-results; -d/--prokfinder-dir is an alias)
    --refdb-dir <dir>                Folder of the reference databases (NCBI, BacDive, GTDB), shared by every
                                     run and downloaded only once (default: \$HOME/prokfinder-refdb, or the
                                     environment variable PROKFINDER_REFDB_DIR)
    -p, --threads <n>                CPUs for seqkit, CheckM2, GTDB-Tk and skDER (default: all, $(nproc 2>/dev/null || echo 4))
    -j, --download-threads <n>       Concurrent genome downloads, 1-30 (default: 10)
    -z, --gzip                       Keep the genome files compressed (.fna.gz)
    --min-contig-length <bp>         Contigs shorter than this are removed (default: 1000)
    --gtdbtk                         Classify the genomes not found in the GTDB tables with GTDB-Tk
                                     (WARNING: needs > 200 GB of RAM); by default they are only listed in
                                     OUTPUT/GENOME-TAXONOMY/NON-CLASSIFIED-GENOMES.txt
    -n, --select-only                Only select the genomes (no download, quality control or taxonomy)
    --pbs                            Run on the cluster: write a PBS job script and submit it with qsub
                                     (resources: PROKFINDER_PBS_QUEUE/PPN/MEM/WALLTIME; default roteamento,
                                     20 CPUs, 20gb; see modules/submit-prokfinder-pbs.sh -h)
    -h, --help                       Show this help message

  Re-running the pipeline:
    Finished work is reused: the reference database is downloaded only once; the selection is reused when the
    filters did not change; downloaded genomes, clean genomes, seqkit/CheckM2/GTDB-Tk results and skDER results
    are kept per genome (or per taxon) and only new genomes are processed. To force a step to run again:
    --force-refdb       Download the NCBI, BacDive and GTDB reference databases again
    --force-selection   Run the genome search again, even if the filters are the same as in the last run
    --force-curation    Discard the clean genomes and every quality/taxonomy/dereplication result
                        (the downloaded genomes are kept)

  Quality thresholds (environment variables): QC_MAX_CONTIGS (1000), QC_MIN_N50 (5000), QC_MAX_N_BASES (100000),
  QC_MIN_COMPLETENESS (75), QC_MAX_CONTAMINATION (5). Required environment variables: CHECKM2DB (and
  GTDBTK_DATA_PATH with --gtdbtk). Optional: NCBI_API_KEY.

  Examples:
    $(basename "$0") -T 45:80 -O rod,coccus -g isolate_genome
    $(basename "$0") -B gram_stain=negative -M anaerobe -C countries.txt -p 32
    $(basename "$0") --gtdb-taxon p__Deinococcota --gtdb-range checkm2_completeness=95: -n
EOF
    exit "${1:-1}"
}

############ READING PARAMETERS ############

QUERY_LINES=()
PROKFINDER_DIR="prokfinder-results"
REFDB_DIR="${PROKFINDER_REFDB_DIR:-${HOME}/prokfinder-refdb}"
ASSEMBLY_SOURCE="genbank"
THREADS="$(nproc 2>/dev/null || echo 4)"
DOWNLOAD_THREADS=10
MIN_CONTIG_LENGTH=1000
GZIP=0
SELECT_ONLY=0
USE_GTDBTK=0
USE_PBS=0
FORCE_REFDB=0
FORCE_SELECTION=0
FORCE_CURATION=0

# Adds one filter to the query ("<KEY><TAB><VALUE>")
add_query() {
    QUERY_LINES+=("$1"$'\t'"$2")
}

# Adds one filter per comma-separated value
add_query_csv() {
    local ITEM
    local -a ITEMS
    IFS=',' read -ra ITEMS <<< "$2"
    for ITEM in "${ITEMS[@]}"; do
        ITEM="$(echo "${ITEM}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
        if [[ -n "${ITEM}" ]]; then
            add_query "$1" "${ITEM}"
        fi
    done
}

# Adds one filter per line of a list file
add_query_file() {
    local LINE
    local N=0
    [[ -f "$2" ]] || die "List file not found: '$2'"
    while IFS= read -r LINE || [[ -n "${LINE}" ]]; do
        LINE="$(echo "${LINE%$'\r'}" | sed -E 's/^[[:space:]]+//; s/[[:space:]]+$//')"
        if [[ -n "${LINE}" ]]; then
            add_query "$1" "${LINE}"
            N=$((N + 1))
        fi
    done < "$2"
    [[ "${N}" -gt 0 ]] || die "No valid entry found in the list file '$2'."
}

# Command line errors: only the error and a hint are printed (the help is shown only with -h/--help)
cli_error() {
    log "error" "$1"
    echo "Use $(basename "$0") -h (or --help) to see the options." >&2
    exit 1
}

# Options without a value, and options that take a value
FLAG_OPTIONS=" -h --help -z --gzip -n --select-only --gtdbtk --pbs --force-refdb --force-selection --force-curation "
VALUE_OPTIONS=" -t --taxid -f --taxid-file -a --accession -A --accession-file --ncbi-taxon -i --isolation-source
 -I --isolation-source-file -c --country -C --country-file -g --genome-type -G --gc --ncbi-range --ncbi-match
 -s --assembly-source -T --temp-opt -P --ph-opt -S --salt-opt --temp-growth --ph-growth --salt-growth
 --bacdive-range -O --cell-shape -M --oxygen -B --bacdive-class --gtdb-taxon --gtdb-range --gtdb-match
 -o --results-dir -d --prokfinder-dir --refdb-dir -p --threads -j --download-threads --min-contig-length "
VALUE_OPTIONS="${VALUE_OPTIONS//$'\n'/ }"

is_value_option() { [[ "${VALUE_OPTIONS}" == *" $1 "* ]]; }
is_known_option() { is_value_option "$1" || [[ "${FLAG_OPTIONS}" == *" $1 "* ]]; }

if [[ $# -eq 0 ]]; then
    cli_error "No arguments provided. At least one selection filter is required."
fi

COMMAND_LINE="$0 $*"
ORIGINAL_ARGS=("$@")

while [[ $# -gt 0 ]]; do
    OPT="$1"

    # "--option=value" is accepted as "--option value"
    if [[ "${OPT}" == --*=* ]]; then
        shift
        set -- "${OPT%%=*}" "${OPT#*=}" "$@"
        OPT="$1"
    fi

    # Options without a value
    case "${OPT}" in
        -h|--help)         usage 0 ;;
        -z|--gzip)         GZIP=1; shift; continue ;;
        -n|--select-only)  SELECT_ONLY=1; shift; continue ;;
        --gtdbtk)          USE_GTDBTK=1; shift; continue ;;
        --force-refdb)     FORCE_REFDB=1; shift; continue ;;
        --force-selection) FORCE_SELECTION=1; shift; continue ;;
        --force-curation)  FORCE_CURATION=1; shift; continue ;;
        --pbs)             USE_PBS=1; shift; continue ;;
    esac

    if ! is_value_option "${OPT}"; then
        if [[ "${OPT}" == -* ]]; then
            cli_error "Unknown option: '${OPT}'."
        fi
        cli_error "Unexpected argument '${OPT}': it is neither an option nor the value of one. Several values of the same option must be separated by commas, without spaces (e.g., -t 562,1280), and a value that contains spaces must be written between double quotes (e.g., -i \"hot spring\")."
    fi

    # The value is missing when the option is the last argument or is followed by another option
    # (a negative number such as "-T -5:10" is a value, not an option)
    if [[ $# -lt 2 ]] || is_known_option "$2"; then
        cli_error "Option ${OPT} requires a value."
    fi
    VAL="$2"
    shift 2

    case "${OPT}" in
        -t|--taxid)                   add_query_csv TAXID "${VAL}" ;;
        -f|--taxid-file)              add_query_file TAXID "${VAL}" ;;
        -a|--accession)               add_query_csv ACCESSION "${VAL}" ;;
        -A|--accession-file)          add_query_file ACCESSION "${VAL}" ;;
        --ncbi-taxon)                 add_query_csv NCBI_TAXON "${VAL}" ;;
        -i|--isolation-source)        add_query_csv ISOLATION_SOURCE "${VAL}" ;;
        -I|--isolation-source-file)   add_query_file ISOLATION_SOURCE "${VAL}" ;;
        -c|--country)                 add_query_csv COUNTRY "${VAL}" ;;
        -C|--country-file)            add_query_file COUNTRY "${VAL}" ;;
        -g|--genome-type)             add_query GENOME_TYPE "${VAL}" ;;
        -G|--gc)                      add_query NCBI_RANGE "ASSEMBLY_STATS_GC_PERCENT=${VAL}" ;;
        --ncbi-range)                 add_query NCBI_RANGE "${VAL}" ;;
        --ncbi-match)                 add_query NCBI_MATCH "${VAL}" ;;
        -s|--assembly-source)         add_query ASSEMBLY_SOURCE "${VAL}"; ASSEMBLY_SOURCE="${VAL,,}" ;;
        -T|--temp-opt)                add_query BACDIVE_RANGE "temperature_opt=${VAL}" ;;
        -P|--ph-opt)                  add_query BACDIVE_RANGE "ph_opt=${VAL}" ;;
        -S|--salt-opt)                add_query BACDIVE_RANGE "halophily_opt=${VAL}" ;;
        --temp-growth)                add_query BACDIVE_RANGE "temperature_preference=${VAL}" ;;
        --ph-growth)                  add_query BACDIVE_RANGE "ph_preference=${VAL}" ;;
        --salt-growth)                add_query BACDIVE_RANGE "halophily_preference=${VAL}" ;;
        --bacdive-range)              add_query BACDIVE_RANGE "${VAL}" ;;
        -O|--cell-shape)              add_query BACDIVE_CLASS "cell_shape=${VAL}" ;;
        -M|--oxygen)                  add_query BACDIVE_CLASS "oxygen_tolerance=${VAL}" ;;
        -B|--bacdive-class)           add_query BACDIVE_CLASS "${VAL}" ;;
        --gtdb-taxon)                 add_query_csv GTDB_TAXON "${VAL}" ;;
        --gtdb-range)                 add_query GTDB_RANGE "${VAL}" ;;
        --gtdb-match)                 add_query GTDB_MATCH "${VAL}" ;;
        -o|--results-dir|-d|--prokfinder-dir) PROKFINDER_DIR="${VAL%/}" ;;
        --refdb-dir)                  REFDB_DIR="${VAL%/}" ;;
        -p|--threads)                 THREADS="${VAL}" ;;
        -j|--download-threads)        DOWNLOAD_THREADS="${VAL}" ;;
        --min-contig-length)          MIN_CONTIG_LENGTH="${VAL}" ;;
    esac
done

[[ "${THREADS}" =~ ^[1-9][0-9]*$ ]] || cli_error "Invalid number of threads '-p ${THREADS}'. It must be a positive integer."
if [[ ! "${DOWNLOAD_THREADS}" =~ ^[0-9]+$ ]] || (( DOWNLOAD_THREADS < 1 || DOWNLOAD_THREADS > 30 )); then
    cli_error "Invalid number of concurrent downloads '-j ${DOWNLOAD_THREADS}'. It must be an integer from 1 to 30."
fi
[[ "${MIN_CONTIG_LENGTH}" =~ ^[0-9]+$ ]] || cli_error "Invalid --min-contig-length '${MIN_CONTIG_LENGTH}'. It must be an integer."

# The query is checked now, before anything is downloaded (the search module has the rules of every filter)
make_temp_dir WORK_DIR
QUERY_TMP="${WORK_DIR}/genome-selection-query.txt"
{
    echo "# Prokfinder genome selection query (${PROKFINDER_VERSION})"
    echo "# Command: ${COMMAND_LINE}"
    echo "# Format: <KEY><TAB><VALUE>; same key in ACCESSION/TAXID/NCBI_TAXON/GTDB_TAXON/ISOLATION_SOURCE/COUNTRY = OR; everything else = AND"
    if [[ "${#QUERY_LINES[@]}" -gt 0 ]]; then
        printf '%s\n' "${QUERY_LINES[@]}"
    fi
} > "${QUERY_TMP}"
bash "${MODULES_DIR}/search-genomes-in-reference-db.sh" -q "${QUERY_TMP}" --validate-only \
    || { echo "Use $(basename "$0") -h (or --help) to see the options." >&2; exit 1; }

############ EXECUTION ON THE CLUSTER (PBS) ############
#   With --pbs the run is not made here: the same options (without --pbs) are written into a PBS job script and
#   submitted with qsub. The job then runs this script on a compute node.

if [[ "${USE_PBS}" -eq 1 ]]; then
    if ! command -v qsub >/dev/null 2>&1; then
        cli_error "Execution on the cluster is not possible: the PBS command 'qsub' was not found. Run Prokfinder without --pbs, or on a machine with PBS."
    fi
    PBS_ARGS=()
    for ARG in "${ORIGINAL_ARGS[@]}"; do
        if [[ "${ARG}" != "--pbs" ]]; then
            PBS_ARGS+=("${ARG}")
        fi
    done
    PBS_STATUS=0
    bash "${MODULES_DIR}/submit-prokfinder-pbs.sh" -- "${PBS_ARGS[@]}" || PBS_STATUS=$?
    exit "${PBS_STATUS}"
fi

############ DIRECTORY STRUCTURE ############

# Reference databases: one folder shared by every run (and every results folder), downloaded only once
REF_DIR="${REFDB_DIR}"
NCBI_DIR="${REF_DIR}/NCBI"
BACDIVE_DIR="${REF_DIR}/BACDIVE"
GTDB_DIR="${REF_DIR}/GTDB"

OUTPUT_DIR="${PROKFINDER_DIR}/OUTPUT"
QUERY_FILE="${OUTPUT_DIR}/genome-selection-query.txt"
METADATA_DIR="${OUTPUT_DIR}/GENOME-METADATA"
METADATA="${METADATA_DIR}/genome-metadata.tsv"
SELECTION_FINGERPRINT="${METADATA_DIR}/.selection-fingerprint"
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
mkdir -p "${REF_DIR}"
log "message" "Results folder: $(absolute_path "${PROKFINDER_DIR}") (log: ${LOG_FILE})"
log "message" "Reference databases folder: $(absolute_path "${REF_DIR}")"

if [[ "${FORCE_CURATION}" -eq 1 ]]; then
    log "message" "Option --force-curation: removing the clean genomes and the quality, taxonomy, hierarchy and final results."
    rm -rf "${CLEAN_DIR}" "${QUALITY_DIR}" "${TAXONOMY_DIR}" "${HIERARCHY_DIR}" "${FINAL_DIR}"
fi

mkdir -p \
    "$(pf_bacdive_brute_dir "${BACDIVE_DIR}")" "$(dirname "$(pf_bacdive_table "${BACDIVE_DIR}")")" \
    "${METADATA_DIR}" "${RAW_DIR}" "${CLEAN_DIR}" "${QUALITY_DIR}" "${TAXONOMY_DIR}" "${HIERARCHY_DIR}" \
    "${FINAL_DIR}/STATISTICS/TABLES" "${FINAL_DIR}/STATISTICS/GRAPHS"
for DOMAIN in "${PF_DOMAINS[@]}"; do
    mkdir -p "$(pf_ncbi_genomes_dir "${NCBI_DIR}" "${DOMAIN}")" "$(pf_ncbi_metadata_dir "${NCBI_DIR}" "${DOMAIN}")" \
             "${GTDB_DIR}/${DOMAIN}-GTDB-TAXONOMY" "${GTDB_DIR}/${DOMAIN}-GTDB-METADATA"
done

############ DEPENDENCIES ############

# Reference databases that must be downloaded in this run (REFDB_NEEDED[<DB>] = 1)
declare -A REFDB_NEEDED=()
check_refdb_needed() {
    local DB DB_DIR
    REFDB_NEEDED=()
    for DB in NCBI BACDIVE GTDB; do
        case "${DB}" in
            NCBI)    DB_DIR="${NCBI_DIR}" ;;
            BACDIVE) DB_DIR="${BACDIVE_DIR}" ;;
            GTDB)    DB_DIR="${GTDB_DIR}" ;;
        esac
        if [[ "${FORCE_REFDB}" -eq 1 || ! -f "$(pf_completed_marker "${DB_DIR}")" ]]; then
            REFDB_NEEDED["${DB}"]=1
        fi
    done
}
check_refdb_needed

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
fi
if [[ -n "${REFDB_NEEDED[GTDB]:-}" ]]; then
    check_program md5sum "coreutils"
fi

if [[ "${SELECT_ONLY}" -eq 0 ]]; then
    check_program datasets "conda install -c conda-forge ncbi-datasets-cli"
    check_program seqkit "conda install -c bioconda seqkit"
    check_program checkm2 "conda install -c bioconda checkm2"
    check_program skder "conda install -c bioconda skder"

    if [[ -z "${CHECKM2DB:-}" ]]; then
        log "error" "CHECKM2DB is not set. It must be the path of the CheckM2 database file, e.g. export CHECKM2DB=/path/CheckM2_database/uniref100.KO.1.dmnd (download it with: checkm2 database --download)."
        DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
    elif [[ ! -f "${CHECKM2DB}" || ! -r "${CHECKM2DB}" ]]; then
        log "error" "CHECKM2DB must point to the database file itself (uniref100.KO.1.dmnd); not found or not readable: ${CHECKM2DB}"
        DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
    fi

    if [[ "${USE_GTDBTK}" -eq 1 ]]; then
        check_program gtdbtk "conda install -c bioconda gtdbtk"
        if [[ -z "${GTDBTK_DATA_PATH:-}" || ! -d "${GTDBTK_DATA_PATH}" ]]; then
            log "error" "GTDBTK_DATA_PATH is not set or is not a folder. It must point to the GTDB-Tk reference package."
            DEPENDENCY_ERRORS=$((DEPENDENCY_ERRORS + 1))
        fi
    fi
fi

if [[ "${DEPENDENCY_ERRORS}" -gt 0 ]]; then
    die "${DEPENDENCY_ERRORS} problem(s) found in the dependencies. Fix them (or use -n to only select the genomes)."
fi

log "message" "Dependencies:"
for TOOL in datasets dataformat taxonkit seqkit checkm2 skder gtdbtk Rscript; do
    if command -v "${TOOL}" >/dev/null 2>&1; then
        case "${TOOL}" in
            seqkit)   log "message" "  seqkit: $(tool_version seqkit version)" ;;
            dataformat) log "message" "  dataformat: $(tool_version dataformat version)" ;;
            taxonkit) log "message" "  taxonkit: $(tool_version taxonkit version)" ;;
            Rscript)  log "message" "  Rscript: $(tool_version Rscript --version)" ;;
            *)        log "message" "  ${TOOL}: $(tool_version "${TOOL}" --version)" ;;
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

############ STEP 0: REFERENCE DATABASE ############

#   The reference folder is shared by every run, so only one run downloads it at a time: a run that finds another
#   one downloading waits for it and then reuses what it downloaded.

FORCE_FLAG=()
if [[ "${FORCE_REFDB}" -eq 1 ]]; then
    FORCE_FLAG=(--force)
    log "message" "Option --force-refdb: the reference databases will be downloaded again."
fi

REFDB_LOCKED=0
if [[ "${#REFDB_NEEDED[@]}" -gt 0 ]]; then
    if command -v flock >/dev/null 2>&1; then
        exec 9> "${REF_DIR}/.download.lock"
        if ! flock -n 9; then
            log "message" "Another Prokfinder run is downloading the reference databases in ${REF_DIR}; waiting for it to finish..."
            flock 9
            check_refdb_needed
        fi
        REFDB_LOCKED=1
    else
        log "warning" "'flock' not found: two Prokfinder runs downloading the reference databases at the same time could conflict."
    fi
fi

if [[ -n "${REFDB_NEEDED[NCBI]:-}" ]]; then
    run_module download-reference-db-ncbi.sh -o "${NCBI_DIR}" --assembly-source all ${FORCE_FLAG[@]+"${FORCE_FLAG[@]}"}
else
    log "message" "NCBI reference data already downloaded ($(head -n 1 "$(pf_completed_marker "${NCBI_DIR}")"))."
fi

if [[ -n "${REFDB_NEEDED[BACDIVE]:-}" ]]; then
    run_module download-reference-db-bacdive.sh -o "${BACDIVE_DIR}" ${FORCE_FLAG[@]+"${FORCE_FLAG[@]}"}
else
    log "message" "BacDive reference data already downloaded ($(head -n 1 "$(pf_completed_marker "${BACDIVE_DIR}")"))."
fi

if [[ -n "${REFDB_NEEDED[GTDB]:-}" ]]; then
    run_module download-reference-db-gtdb.sh -o "${GTDB_DIR}" ${FORCE_FLAG[@]+"${FORCE_FLAG[@]}"}
else
    log "message" "GTDB reference data already downloaded ($(head -n 1 "$(pf_completed_marker "${GTDB_DIR}")"))."
fi

if [[ "${REFDB_LOCKED}" -eq 1 ]]; then
    flock -u 9
    exec 9>&-
fi

############ STEP 1: GENOME SEARCH ############
#   The selection is reused when the query and the reference databases are the same as in the last run.

CURRENT_FINGERPRINT="$(
    grep -v '^#' "${QUERY_TMP}"
    for DB_DIR in "${NCBI_DIR}" "${BACDIVE_DIR}" "${GTDB_DIR}"; do
        echo "${DB_DIR}: $(file_mtime "$(pf_completed_marker "${DB_DIR}")")"
    done
    echo "BACDIVE_TABLE: $(file_mtime "$(pf_bacdive_table "${BACDIVE_DIR}")")"
)"

if [[ "${FORCE_SELECTION}" -eq 0 && -s "${METADATA}" && -f "${SELECTION_FINGERPRINT}" \
      && "$(cat "${SELECTION_FINGERPRINT}")" == "${CURRENT_FINGERPRINT}" ]]; then
    log "message" "The selection filters and reference databases did not change since the last run: reusing ${METADATA} (use --force-selection to search again)."
else
    if [[ "${FORCE_SELECTION}" -eq 1 ]]; then
        log "message" "Option --force-selection: running the genome search again."
    fi
    log "message" "Selection query:"
    while IFS= read -r LINE; do
        log "message" "    ${LINE//$'\t'/ = }"
    done < <(grep -v '^#' "${QUERY_TMP}")

    rm -f "${SELECTION_FINGERPRINT}"
    run_module search-genomes-in-reference-db.sh -q "${QUERY_TMP}" -r "${REF_DIR}" -o "${METADATA}"
    cp "${QUERY_TMP}" "${QUERY_FILE}"
    printf '%s\n' "${CURRENT_FINGERPRINT}" > "${SELECTION_FINGERPRINT}"
fi

SELECTED_LIST="${METADATA_DIR}/selected-accessions.txt"
tail -n +2 "${METADATA}" | cut -f1 > "${SELECTED_LIST}"
N_SELECTED="$(count_lines "${SELECTED_LIST}")"
log "message" "Selected genome assemblies: ${N_SELECTED} (metadata: ${METADATA})"

if [[ "${N_SELECTED}" -eq 0 ]]; then
    log "message" "No genome assembly matches all the criteria, so there is nothing to download. The counts above show the filter that removed the last genomes."
    exit 0
fi

if [[ "${SELECT_ONLY}" -eq 1 ]]; then
    log "message" "Option -n: only the selection was made (no download, quality control or taxonomy)."
    exit 0
fi

GZIP_FLAG=()
if [[ "${GZIP}" -eq 1 ]]; then
    GZIP_FLAG=(-z)
fi

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

############ STEP 2: GENOME DOWNLOAD ############

run_module download-selected-genomes-from-db.sh -l "${SELECTED_LIST}" -r "${NCBI_DIR}" -o "${RAW_DIR}" \
    -t "${GENOME_FILES_DIR}/genome-download-status.tsv" -j "${DOWNLOAD_THREADS}" ${GZIP_FLAG[@]+"${GZIP_FLAG[@]}"}
add_metadata_columns "${METADATA}" "${GENOME_FILES_DIR}/genome-download-status.tsv"

DOWNLOADED_LIST="${WORK_DIR}/downloaded.txt"
accessions_where "${DOWNLOADED_LIST}" '$DOWNLOAD_STATUS == "DOWNLOADED"' DOWNLOAD_STATUS
if [[ ! -s "${DOWNLOADED_LIST}" ]]; then
    finish_early "No selected genome could be downloaded."
fi

############ STEP 3: CONTIG FILTERING ############

run_module filter-low-quality-contigs-from-genomes.sh -l "${DOWNLOADED_LIST}" -i "${RAW_DIR}" -o "${CLEAN_DIR}" \
    -t "${QUALITY_DIR}/contig-filtering.tsv" -m "${MIN_CONTIG_LENGTH}" -p "${THREADS}" ${GZIP_FLAG[@]+"${GZIP_FLAG[@]}"}
add_metadata_columns "${METADATA}" "${QUALITY_DIR}/contig-filtering.tsv"

CLEAN_LIST="${WORK_DIR}/clean.txt"
accessions_where "${CLEAN_LIST}" '$CONTIG_FILTER_STATUS == "OK"' CONTIG_FILTER_STATUS
if [[ ! -s "${CLEAN_LIST}" ]]; then
    finish_early "No genome has contigs left after the contig filtering."
fi

############ STEP 4: QUALITY CONTROL ############

run_module evaluate-and-filter-genome-by-qc.sh -l "${CLEAN_LIST}" -i "${CLEAN_DIR}" -o "${QUALITY_DIR}" -p "${THREADS}"
add_metadata_columns "${METADATA}" "${QUALITY_DIR}/genome-quality.tsv"

HIGH_QUALITY_LIST="${WORK_DIR}/high_quality.txt"
accessions_where "${HIGH_QUALITY_LIST}" '$QC_STATUS == "PASSED"' QC_STATUS
if [[ ! -s "${HIGH_QUALITY_LIST}" ]]; then
    finish_early "No genome passed the quality filters."
fi

############ STEP 5: GTDB TAXONOMY ############

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

############ STEP 6: TAXONOMIC DIRECTORY HIERARCHY ############

run_module build-taxonomic-dir-hierarchy.sh -m "${METADATA}" -l "${CLASSIFIED_LIST}" -o "${HIERARCHY_DIR}"
add_metadata_columns "${METADATA}" "${HIERARCHY_DIR}/genome-dir-hierarchy.tsv"

############ STEP 7: DEREPLICATION ############

run_module dereplicate-genomes-skder.sh -d "${HIERARCHY_DIR}" -m "${METADATA}" \
    -t "${HIERARCHY_DIR}/genome-dereplication.tsv" -p "${THREADS}"
add_metadata_columns "${METADATA}" "${HIERARCHY_DIR}/genome-dereplication.tsv"

############ STEP 8: FINAL RESULTS ############

run_module summarize-final-results.sh -m "${METADATA}" -r "${NCBI_DIR}" -s "${ASSEMBLY_SOURCE}" -o "${FINAL_DIR}"

log "message" "Prokfinder finished successfully. Genome metadata (every step): ${METADATA}"
log "message" "High-quality representative genomes: ${FINAL_DIR}/high-quality-selected-genomes.txt"
