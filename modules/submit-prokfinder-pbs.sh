#!/usr/bin/env bash
#
# submit-prokfinder-pbs.sh - Prokfinder module: execution on a cluster with the PBS queue manager
#
# Writes a PBS submission script (.pbs) that loads the environment modules of the cluster and runs
# main-prokfinder.sh with the given options, then sends it to the queue with qsub. It is called by
# "main-prokfinder.sh --pbs ..." and can also be run directly.
#
# Usage:
#   submit-prokfinder-pbs.sh [--dry-run] -- <main-prokfinder.sh options>   submit a Prokfinder run
#   submit-prokfinder-pbs.sh --test [--dry-run]                             submit a short job that only checks
#                                                                           the modules and the programs
#
# The .pbs script and the job log are written to <prokfinder_dir>/PBS-JOBS/ (<prokfinder_dir> = the -d option
# of Prokfinder, -o/--results-dir, default ./prokfinder-results). The job runs in the folder from where it was submitted.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROKFINDER_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

# ---------- Defaults (each one can be changed through the environment variable of the same name) ----------
PROKFINDER_PBS_QUEUE="${PROKFINDER_PBS_QUEUE:-roteamento}"
PROKFINDER_PBS_PPN="${PROKFINDER_PBS_PPN:-20}"
PROKFINDER_PBS_MEM="${PROKFINDER_PBS_MEM:-20gb}"
PROKFINDER_PBS_JOB_NAME="${PROKFINDER_PBS_JOB_NAME:-prokfinder}"
PROKFINDER_PBS_MODULES="${PROKFINDER_PBS_MODULES:-seqkit checkm2 ncbi_datasets gtdbtk skder taxonkit}"
PROKFINDER_PBS_OPTIONAL_MODULES="${PROKFINDER_PBS_OPTIONAL_MODULES:-r}"

usage() {
    cat << EOF
Usage:
  $(basename "$0") [--dry-run] -- <main-prokfinder.sh options>
  $(basename "$0") --test [--dry-run]
  (usually called through: main-prokfinder.sh --pbs <options>)

Writes a PBS submission script that loads the cluster modules and runs Prokfinder, and submits it with qsub.

Options:
  --test       Submit a short job that only loads the modules and prints the version/help of every program
               (seqkit, checkm2, datasets, dataformat, gtdbtk, skder, taxonkit, R). Default resources of the test:
               queue qtime, 1 CPU, 2gb, 00:15:00
  --dry-run    Only write the .pbs script (do not submit it)
  -h, --help   Show this help message

Resources (environment variables; the first four are also options of main-prokfinder.sh):
  PROKFINDER_PBS_QUEUE      --pbs-queue: roteamento, qtime or high_mem (default: roteamento)
  PROKFINDER_PBS_PPN        --pbs-cpus: CPUs of the job (default: 20). Prokfinder uses all of them unless -p is given
  PROKFINDER_PBS_MEM        --pbs-mem: RAM of the job, e.g. 20gb, 20000mb or 20 (GB) (default: 20gb)
  PROKFINDER_PBS_WALLTIME   --pbs-walltime: maximum run time HH:MM:SS (default: 2160:00:00; 72:00:00 for qtime)
  PROKFINDER_PBS_JOB_NAME   Job name (default: prokfinder)
  PROKFINDER_PBS_MODULES    Modules loaded in the job (default: "seqkit checkm2 ncbi_datasets gtdbtk skder taxonkit")
  PROKFINDER_PBS_OPTIONAL_MODULES  Modules loaded if available (default: "r", for the graphs)

Queues of the UFV cluster (one node per job):
  roteamento  routes to "work": up to 24024 MB of RAM per job; at most 96 CPUs and 96 GB per user;
              walltime up to 2160:00:00 (3 months)
  qtime       short jobs, routes to "workt": up to 24024 MB per job; walltime REQUIRED, at most 72:00:00
  high_mem    at least 25 GB of RAM per job; at most 64 CPUs and 512 GB per user; 2 running jobs per user;
              walltime up to 2160:00:00. Needed for GTDB-Tk (--gtdbtk needs > 200 GB of RAM)

Examples:
  main-prokfinder.sh --pbs -T 45:80 -O rod
  main-prokfinder.sh --pbs --pbs-cpus 32 --pbs-mem 24000mb -T 45:80
  main-prokfinder.sh --pbs --pbs-queue high_mem --pbs-mem 256gb --pbs-cpus 32 -T 45:80 --gtdbtk
  PROKFINDER_PBS_QUEUE=qtime PROKFINDER_PBS_WALLTIME=24:00:00 main-prokfinder.sh --pbs -T 45:80 -n
  $(basename "$0") --test
EOF
    exit "${1:-0}"
}

pbs_error() {
    log "error" "$1"
    echo "Use $(basename "$0") -h to see the options." >&2
    exit 1
}

TEST=0
DRY_RUN=0
MAIN_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --test)    TEST=1; shift ;;
        --dry-run) DRY_RUN=1; shift ;;
        -h|--help) usage 0 ;;
        --)        shift; MAIN_ARGS=("$@"); break ;;
        *)         pbs_error "Unknown option: '$1'. The Prokfinder options must come after '--'." ;;
    esac
done

if [[ "${TEST}" -eq 0 && "${#MAIN_ARGS[@]}" -eq 0 ]]; then
    pbs_error "No Prokfinder options given (they must come after '--'), and --test was not used."
fi

if [[ "${TEST}" -eq 1 ]]; then
    PROKFINDER_PBS_QUEUE="${PROKFINDER_PBS_TEST_QUEUE:-qtime}"
    PROKFINDER_PBS_PPN=1
    PROKFINDER_PBS_MEM=2gb
    PROKFINDER_PBS_WALLTIME="00:15:00"
    PROKFINDER_PBS_JOB_NAME="prokfinder_test"
fi

if [[ -z "${PROKFINDER_PBS_WALLTIME:-}" ]]; then
    if [[ "${PROKFINDER_PBS_QUEUE}" == "qtime" ]]; then
        PROKFINDER_PBS_WALLTIME="72:00:00"
    else
        PROKFINDER_PBS_WALLTIME="2160:00:00"
    fi
fi

############ CHECKING THE RESOURCES ############

# Prints the amount of memory in MB ("20gb" or "20" -> 20480, "24000mb" -> 24000); returns 1 if the format is invalid
mem_to_mb() {
    local VALUE="${1,,}"
    if [[ "${VALUE}" =~ ^([0-9]+)(gb|g)?$ ]]; then
        echo $(( BASH_REMATCH[1] * 1024 ))
    elif [[ "${VALUE}" =~ ^([0-9]+)(mb|m)$ ]]; then
        echo "${BASH_REMATCH[1]}"
    else
        return 1
    fi
}

# Prints a walltime HH:MM:SS in hours (rounded up); returns 1 if the format is invalid
walltime_hours() {
    if [[ "$1" =~ ^([0-9]+):([0-5][0-9]):([0-5][0-9])$ ]]; then
        local SECONDS_TOTAL=$(( 10#${BASH_REMATCH[1]} * 3600 + 10#${BASH_REMATCH[2]} * 60 + 10#${BASH_REMATCH[3]} ))
        echo $(( (SECONDS_TOTAL + 3599) / 3600 ))
    else
        return 1
    fi
}

[[ "${PROKFINDER_PBS_PPN}" =~ ^[1-9][0-9]*$ ]] || pbs_error "The CPUs of the job (--pbs-cpus / PROKFINDER_PBS_PPN) must be a positive integer (got '${PROKFINDER_PBS_PPN}')."
MEM_MB="$(mem_to_mb "${PROKFINDER_PBS_MEM}")" || pbs_error "The RAM of the job (--pbs-mem / PROKFINDER_PBS_MEM) must be like 20gb, 20000mb or 20 (GB) (got '${PROKFINDER_PBS_MEM}')."
# PBS needs the unit: a number alone is in GB
[[ "${PROKFINDER_PBS_MEM}" =~ ^[0-9]+$ ]] && PROKFINDER_PBS_MEM="${PROKFINDER_PBS_MEM}gb"
HOURS="$(walltime_hours "${PROKFINDER_PBS_WALLTIME}")" || pbs_error "The walltime of the job (--pbs-walltime / PROKFINDER_PBS_WALLTIME) must be HH:MM:SS (got '${PROKFINDER_PBS_WALLTIME}')."
[[ "${PROKFINDER_PBS_JOB_NAME}" =~ ^[A-Za-z][A-Za-z0-9_.-]{0,14}$ ]] \
    || pbs_error "PROKFINDER_PBS_JOB_NAME must start with a letter and have at most 15 letters, digits, '_', '.' or '-' (got '${PROKFINDER_PBS_JOB_NAME}')."

case "${PROKFINDER_PBS_QUEUE}" in
    roteamento|qtime)
        (( MEM_MB <= 24024 )) || pbs_error "Queue ${PROKFINDER_PBS_QUEUE} accepts at most 24024 MB of RAM per job (asked: ${PROKFINDER_PBS_MEM}). Use --pbs-queue high_mem for jobs with 25 GB or more."
        (( PROKFINDER_PBS_PPN <= 96 )) || pbs_error "Queue ${PROKFINDER_PBS_QUEUE} accepts at most 96 CPUs per user (asked: ${PROKFINDER_PBS_PPN})."
        if [[ "${PROKFINDER_PBS_QUEUE}" == "qtime" ]]; then
            (( HOURS <= 72 )) || pbs_error "Queue qtime accepts a walltime of at most 72:00:00 (asked: ${PROKFINDER_PBS_WALLTIME})."
        else
            (( HOURS <= 2160 )) || pbs_error "Queue roteamento accepts a walltime of at most 2160:00:00 (asked: ${PROKFINDER_PBS_WALLTIME})."
        fi
        ;;
    high_mem)
        (( MEM_MB >= 25600 )) || pbs_error "Queue high_mem is for jobs with at least 25 GB of RAM (asked: ${PROKFINDER_PBS_MEM})."
        (( MEM_MB <= 524288 )) || pbs_error "Queue high_mem accepts at most 512 GB of RAM per user (asked: ${PROKFINDER_PBS_MEM})."
        (( PROKFINDER_PBS_PPN <= 64 )) || pbs_error "Queue high_mem accepts at most 64 CPUs per user (asked: ${PROKFINDER_PBS_PPN})."
        (( HOURS <= 2160 )) || pbs_error "Queue high_mem accepts a walltime of at most 2160:00:00 (asked: ${PROKFINDER_PBS_WALLTIME})."
        ;;
    *)
        log "warning" "Queue '${PROKFINDER_PBS_QUEUE}' is not one of the known queues (roteamento, qtime, high_mem); its limits were not checked."
        ;;
esac

if [[ "${DRY_RUN}" -eq 0 ]] && ! command -v qsub >/dev/null 2>&1; then
    pbs_error "Execution on the cluster is not possible: the PBS command 'qsub' was not found. Run Prokfinder without --pbs, or on a machine with PBS (use --dry-run to only write the .pbs script)."
fi

############ OPTIONS OF THE PROKFINDER RUN ############

PROKFINDER_DIR="prokfinder-results"
USER_THREADS=""
USES_GTDBTK=0
for (( i = 0; i < ${#MAIN_ARGS[@]}; i++ )); do
    ARG="${MAIN_ARGS[$i]}"
    NEXT="${MAIN_ARGS[$((i + 1))]:-}"
    case "${ARG}" in
        -o|--results-dir|-d|--prokfinder-dir) PROKFINDER_DIR="${NEXT}" ;;
        --results-dir=*|--prokfinder-dir=*) PROKFINDER_DIR="${ARG#*=}" ;;
        -p|--threads)        USER_THREADS="${NEXT}" ;;
        --threads=*)         USER_THREADS="${ARG#*=}" ;;
        --gtdbtk)            USES_GTDBTK=1 ;;
    esac
done

# Prokfinder uses the CPUs of the job, unless the user chose another number
if [[ "${TEST}" -eq 0 ]]; then
    if [[ -z "${USER_THREADS}" ]]; then
        MAIN_ARGS+=(-p "${PROKFINDER_PBS_PPN}")
    elif [[ "${USER_THREADS}" =~ ^[0-9]+$ ]] && (( USER_THREADS > PROKFINDER_PBS_PPN )); then
        log "warning" "-p ${USER_THREADS} is larger than the CPUs of the job (--pbs-cpus ${PROKFINDER_PBS_PPN}); the programs will compete for the CPUs."
    fi
    if [[ "${USES_GTDBTK}" -eq 1 && ( "${PROKFINDER_PBS_QUEUE}" != "high_mem" || MEM_MB -lt 204800 ) ]]; then
        log "warning" "--gtdbtk needs more than 200 GB of RAM, but the job asks for ${PROKFINDER_PBS_MEM} in queue ${PROKFINDER_PBS_QUEUE}. Consider: --pbs-queue high_mem --pbs-mem 256gb"
    fi
fi

mkdir -p "${PROKFINDER_DIR}/PBS-JOBS"
JOBS_DIR="$(absolute_path "${PROKFINDER_DIR}/PBS-JOBS")"
STAMP="$(date +%Y%m%d_%H%M%S)"
PBS_FILE="${JOBS_DIR}/${PROKFINDER_PBS_JOB_NAME}_${STAMP}.pbs"
JOB_LOG="${JOBS_DIR}/${PROKFINDER_PBS_JOB_NAME}_${STAMP}.log"

# Prokfinder settings of the current environment that are copied into the job
# (NCBI_API_KEY is passed by qsub at submission and never written to the .pbs file)
EXPORT_VARS=(
    QC_MAX_CONTIGS QC_MIN_N50 QC_MAX_N_BASES QC_MIN_COMPLETENESS QC_MAX_CONTAMINATION CHECKM2_BATCH_SIZE
    SKDER_MODE SKDER_PERCENT_IDENTITY SKDER_ALIGNED_FRACTION SKDER_EXTRA_ARGS
    GTDBTK_EXTRA_ARGS NCBI_FIRST_YEAR NCBI_MAX_ATTEMPTS NCBI_RETRY_WAIT REHYDRATE_MAX_ATTEMPTS REHYDRATE_RETRY_WAIT
    CHECKM2DB GTDBTK_DATA_PATH PROKFINDER_REFDB_DIR
)

############ WRITING THE PBS SCRIPT ############

{
    cat << EOF
#!/usr/bin/env bash

#PBS -N ${PROKFINDER_PBS_JOB_NAME}
#PBS -l nodes=1:ppn=${PROKFINDER_PBS_PPN},mem=${PROKFINDER_PBS_MEM},walltime=${PROKFINDER_PBS_WALLTIME}
#PBS -q ${PROKFINDER_PBS_QUEUE}
#PBS -j oe
#PBS -o ${JOB_LOG}

# Written by Prokfinder (modules/submit-prokfinder-pbs.sh) on $(date +"%Y-%m-%d %H:%M:%S")
# Resources: change them with the main-prokfinder.sh options --pbs-queue, --pbs-cpus, --pbs-mem and --pbs-walltime

# Change to the folder from where the job was submitted
cd "\${PBS_O_WORKDIR}" || exit 1

# Enable and load the environment modules
source /etc/profile.d/modules.sh
module load ${PROKFINDER_PBS_MODULES} || { echo "ERROR: could not load the modules: ${PROKFINDER_PBS_MODULES}" >&2; exit 1; }
EOF
    for MODULE in ${PROKFINDER_PBS_OPTIONAL_MODULES}; do
        echo "module load ${MODULE} 2>/dev/null || echo \"WARNING: optional module '${MODULE}' could not be loaded\" >&2"
    done
    echo

    if [[ "${TEST}" -eq 0 ]]; then
        echo "# Prokfinder settings of the submission environment"
        for VAR in "${EXPORT_VARS[@]}"; do
            if [[ -n "${!VAR:-}" ]]; then
                printf 'export %s=%q\n' "${VAR}" "${!VAR}"
            fi
        done
        echo
    fi

    cat << 'EOF'
echo "Job ${PBS_JOBID:-?} (${PBS_JOBNAME:-?}) on $(hostname) - started $(date +"%Y-%m-%d %H:%M:%S")"
module list 2>&1
echo
EOF

    if [[ "${TEST}" -eq 1 ]]; then
        cat << 'EOF'
# Test: version and help of every program used by Prokfinder
STATUS=0
check() {
    local NAME="$1"
    shift
    echo "==================== ${NAME}: $* ===================="
    if ! command -v "${NAME}" >/dev/null 2>&1; then
        echo "MISSING: ${NAME} is not in the PATH"
        STATUS=1
        return 0
    fi
    "$@" 2>&1 | head -n 60
    echo
}
check seqkit seqkit version
check checkm2 checkm2 --version
check checkm2 checkm2 predict --help
check datasets datasets --version
check dataformat dataformat version
check taxonkit taxonkit version
check gtdbtk gtdbtk --version
check gtdbtk gtdbtk classify_wf --help
check skder skder --version
check skder skder --help
check Rscript Rscript --version
for TOOL in jq gawk unzip curl md5sum; do
    printf '%-10s %s\n' "${TOOL}" "$(command -v "${TOOL}" || echo MISSING)"
done
echo "CHECKM2DB=${CHECKM2DB:-unset}"
echo "GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH:-unset}"
echo "bash ${BASH_VERSION}"
echo "Finished $(date +"%Y-%m-%d %H:%M:%S") - status ${STATUS}"
exit "${STATUS}"
EOF
    else
        echo "# Prokfinder run"
        printf 'bash %q' "${PROKFINDER_ROOT}/main-prokfinder.sh"
        printf ' %q' "${MAIN_ARGS[@]}"
        echo
        cat << 'EOF'
STATUS=$?
echo "Finished $(date +"%Y-%m-%d %H:%M:%S") - exit status ${STATUS}"
exit "${STATUS}"
EOF
    fi
} > "${PBS_FILE}"

log "message" "PBS script written: ${PBS_FILE}"
log "message" "Resources: queue ${PROKFINDER_PBS_QUEUE}, ${PROKFINDER_PBS_PPN} CPUs, ${PROKFINDER_PBS_MEM}, walltime ${PROKFINDER_PBS_WALLTIME}"

if [[ "${DRY_RUN}" -eq 1 ]]; then
    log "message" "Option --dry-run: the job was not submitted. Submit it with: qsub ${PBS_FILE}"
    exit 0
fi

QSUB_ARGS=()
if [[ -n "${NCBI_API_KEY:-}" ]]; then
    QSUB_ARGS=(-v NCBI_API_KEY)
fi

JOB_ID="$(qsub ${QSUB_ARGS[@]+"${QSUB_ARGS[@]}"} "${PBS_FILE}")" || pbs_error "qsub refused the job (see the message above). The script was kept: ${PBS_FILE}"

log "message" "Job submitted: ${JOB_ID}"
log "message" "Follow it with: qstat -u ${USER:-$(whoami)}    |    log of the job: ${JOB_LOG}"
