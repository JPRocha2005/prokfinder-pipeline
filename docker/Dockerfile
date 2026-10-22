# Prokfinder container image
#
# Build (from the root of the repository):
#     docker build -f docker/Dockerfile -t prokfinder .
#     docker build -f docker/Dockerfile --build-arg WITH_GTDBTK=0 -t prokfinder .   # without GTDB-Tk (smaller)
# Run: docker/prokfinder-docker.sh <sub-command> [options]   (see the README, "Instalação via Docker")
#
# The programs are installed with micromamba (bioconda/conda-forge). CheckM2, BUSCO and GTDB-Tk have dependencies
# that conflict with each other, so each one has its own environment; a small wrapper in /usr/local/bin runs each of
# them from its environment, so Prokfinder calls them as usual. The reference databases (Prokfinder, CheckM2,
# GTDB-Tk, BUSCO lineages) are not in the image: they are mounted from the host.

FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends gawk jq curl unzip ca-certificates procps \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER docker/env-*.yml /tmp/envs/
RUN micromamba install -y -n base -f /tmp/envs/env-prokfinder.yml \
    && micromamba create -y -n checkm2 -f /tmp/envs/env-checkm2.yml \
    && micromamba create -y -n busco -f /tmp/envs/env-busco.yml \
    && micromamba clean -a -y

# GTDB-Tk is optional (ANI of the genomes that are not in the GTDB tables, and the tree with download-wf
# --use-gtdbtk-tree; its reference package is mounted at run time)
ARG WITH_GTDBTK=1
RUN if [ "${WITH_GTDBTK}" = "1" ]; then \
        micromamba create -y -n gtdbtk -f /tmp/envs/env-gtdbtk.yml && micromamba clean -a -y; \
    fi

USER root
# Wrappers: run each program from its own environment (the environment name is also the program name; the variable
# is not called ENV_NAME, which the micromamba image defines and Docker would replace)
RUN set -e; for TOOL_ENV in checkm2 busco gtdbtk; do \
        if [ -d "/opt/conda/envs/${TOOL_ENV}" ]; then \
            printf '#!/usr/bin/env bash\nexport PATH="/opt/conda/envs/%s/bin:${PATH}"\nexec %s "$@"\n' \
                "${TOOL_ENV}" "${TOOL_ENV}" > "/usr/local/bin/${TOOL_ENV}"; \
            chmod 755 "/usr/local/bin/${TOOL_ENV}"; \
        fi; \
    done

COPY . /opt/prokfinder
RUN chmod 755 /opt/prokfinder/main-prokfinder.sh /opt/prokfinder/modules/*.sh /opt/prokfinder/modules/*.py \
        /opt/prokfinder/modules/plot-statistics.R \
    && mkdir -p /home/prokfinder/prokfinder-refdb && chmod -R 777 /home/prokfinder

# HOME holds the default reference databases folder ($HOME/prokfinder-refdb), mounted from the host
ENV HOME=/home/prokfinder \
    PATH=/opt/prokfinder:/opt/conda/bin:${PATH}
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "main-prokfinder.sh"]
CMD ["--help"]
