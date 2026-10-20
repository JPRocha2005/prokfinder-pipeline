#!/usr/bin/env bash
#
# prokfinder-docker.sh - runs Prokfinder from its Docker image, with the same command line as main-prokfinder.sh:
#
#     docker/prokfinder-docker.sh setup
#     docker/prokfinder-docker.sh search-genomes -T 45:80 AND -O rod download-wf --min-completeness 90
#
# Folders mounted in the container (at the same paths, so every path of the command line keeps working):
#   - the current folder (results: -o, list files of the filters...)
#   - the reference databases folder: --refdb-dir if given, otherwise $HOME/prokfinder-refdb
#   - the results folder of -o/--results-dir, when it is outside the current folder
#   - the folders of CHECKM2DB and GTDBTK_DATA_PATH, when they are set (they are passed to the container), and
#     NCBI_API_KEY is passed when it is set
# The container runs with the user and group of the caller, so the results belong to the caller.
#
# The image is "prokfinder" (build it with: docker build -f docker/Dockerfile -t prokfinder .); another one can be
# chosen with --image <name> as the first argument.

set -euo pipefail

IMAGE="prokfinder"
if [[ "${1:-}" == "--image" ]]; then
    IMAGE="$2"
    shift 2
fi

command -v docker >/dev/null 2>&1 || { echo "ERROR: docker was not found in the PATH." >&2; exit 1; }
docker image inspect "${IMAGE}" >/dev/null 2>&1 \
    || { echo "ERROR: Docker image '${IMAGE}' not found. Build it with: docker build -f docker/Dockerfile -t ${IMAGE} ." >&2; exit 1; }

# Prints the absolute path of a file or folder (it does not need to exist)
abs_path() {
    case "$1" in
        /*) echo "$1" ;;
        *)  echo "${PWD}/$1" ;;
    esac
}

REFDB=""
RESULTS=""
ARGS=("$@")
for (( I = 0; I < ${#ARGS[@]}; I++ )); do
    case "${ARGS[$I]}" in
        --refdb-dir)        REFDB="${ARGS[$((I + 1))]:-}" ;;
        --refdb-dir=*)      REFDB="${ARGS[$I]#*=}" ;;
        -o|--results-dir)   RESULTS="${ARGS[$((I + 1))]:-}" ;;
        --results-dir=*)    RESULTS="${ARGS[$I]#*=}" ;;
    esac
done

MOUNTS=(-v "${PWD}:${PWD}" -w "${PWD}")
ENVS=()

# Reference databases: the folder given with --refdb-dir (same path), or $HOME/prokfinder-refdb as the default folder
# of the container (/home/prokfinder/prokfinder-refdb)
if [[ -n "${REFDB}" ]]; then
    REFDB="$(abs_path "${REFDB}")"
    mkdir -p "${REFDB}"
    MOUNTS+=(-v "${REFDB}:${REFDB}")
else
    mkdir -p "${HOME}/prokfinder-refdb"
    MOUNTS+=(-v "${HOME}/prokfinder-refdb:/home/prokfinder/prokfinder-refdb")
fi

if [[ -n "${RESULTS}" ]]; then
    RESULTS="$(abs_path "${RESULTS}")"
    if [[ "${RESULTS}" != "${PWD}"/* ]]; then
        mkdir -p "${RESULTS}"
        MOUNTS+=(-v "${RESULTS}:${RESULTS}")
    fi
fi

if [[ -n "${CHECKM2DB:-}" ]]; then
    CHECKM2DB="$(abs_path "${CHECKM2DB}")"
    MOUNTS+=(-v "$(dirname "${CHECKM2DB}"):$(dirname "${CHECKM2DB}"):ro")
    ENVS+=(-e "CHECKM2DB=${CHECKM2DB}")
fi
if [[ -n "${GTDBTK_DATA_PATH:-}" ]]; then
    GTDBTK_DATA_PATH="$(abs_path "${GTDBTK_DATA_PATH}")"
    MOUNTS+=(-v "${GTDBTK_DATA_PATH}:${GTDBTK_DATA_PATH}:ro")
    ENVS+=(-e "GTDBTK_DATA_PATH=${GTDBTK_DATA_PATH}")
fi
[[ -z "${NCBI_API_KEY:-}" ]] || ENVS+=(-e NCBI_API_KEY)

TTY=()
[[ -t 0 && -t 1 ]] && TTY=(-it)

exec docker run --rm ${TTY[@]+"${TTY[@]}"} --user "$(id -u):$(id -g)" \
    "${MOUNTS[@]}" ${ENVS[@]+"${ENVS[@]}"} "${IMAGE}" "$@"
