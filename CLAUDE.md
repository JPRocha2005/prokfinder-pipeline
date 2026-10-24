# Prokfinder: guide for AI coding agents

Prokfinder is a Bash pipeline (with one Python and one R helper) that selects, downloads and curates prokaryote
genomes (Bacteria and Archaea) from NCBI, combining NCBI assembly metadata, BacDive phenotypes and the GTDB
taxonomy. It is developed at LBCM-UFV (Universidade Federal de Viçosa) to build a non-redundant database of
high-quality thermophile genomes. This file is the starting point for any agent working on the code; the user
documentation is `README.md` (in Portuguese).

Index files kept in the repository:

| File | What it is | How it is kept up to date |
|---|---|---|
| `CLAUDE.md` (this file) | architecture, conventions, how to test | by hand, when a convention or the architecture changes |
| `catalog.json` | every sub-command, option, module, shared function, dependency and table format | **generated**: `python3 tools/build-catalog.py` (never edit by hand) |
| `CHANGELOG.md` | user-visible changes per version (`Unreleased` until a release) | by hand, in every change that a user would notice |
| `llms.txt` | short overview and index of the documentation | by hand, rarely |
| `.claude/skills/update-docs/SKILL.md` | the checklist to run after changing the code | - |

## Architecture

```
main-prokfinder.sh          controller: reads the command line, checks the dependencies and the reference
                            databases, then calls the modules in order. Sub-commands (run in this order):
                              setup           reference databases (NCBI, BacDive, GTDB)
                              search-genomes  genome selection (filters joined with AND/OR)
                              download-wf     download -> curation -> taxonomy -> dereplication -> summary
                              run-on-cluster  writes a PBS job that runs the other sub-commands
lib/common.sh               sourced by every script: log/die, temporary folders, TSV helpers (format_tsv,
                            add_metadata_columns), per-genome caches, the NULL rule (PF_AWK_NULL), paths of the
                            reference databases (pf_*), table formats (PF_REFDB_FORMAT_*), run_in_own_env
modules/                    one script per step; each one runs alone too (-h shows its usage)
  download-reference-db-{ncbi,bacdive,gtdb}.sh   setup (0.1-0.3)
  search-genomes-in-reference-db.sh              1  selection (reads genome-selection-query.txt)
  download-selected-genomes-from-db.sh           2  datasets rehydrate
  filter-low-quality-contigs-from-genomes.sh     3  seqkit
  evaluate-genome-structure-seqkit.sh            4  contigs, N50, N bases
  predict-genes-pyrodigal.sh (+ .py)             5  genes predicted once, used by every later step
  evaluate-genome-completeness.sh                6  CheckM2 --genes (filter) + BUSCO proteins (reported)
  assign-taxonomy-gtdb.sh                        7  GTDB table -> GTDB-Tk ani_rep -> tree (--use-gtdbtk-tree)
  build-taxonomic-dir-hierarchy.sh               8  Domain/.../Species folders (lists of paths)
  dereplicate-genomes.sh                         9  skDER or CiDDER inside each leaf
  summarize-final-results.sh (+ plot-statistics.R) 10 tables and graphs
  submit-prokfinder-pbs.sh                       run-on-cluster
docker/                     Dockerfile (micromamba; CheckM2, BUSCO and GTDB-Tk in their own environments) and
                            prokfinder-docker.sh; image published by .github/workflows/docker-publish.yml
tools/build-catalog.py      writes catalog.json
```

Data flow: every step adds columns to one table, `OUTPUT/GENOME-METADATA/genome-metadata.tsv` (one row per genome,
first column ASSEMBLY_ACCESSION), with `add_metadata_columns`. Each step reads the accessions that passed the
previous one (`accessions_where` in main-prokfinder.sh).

Reference databases (`~/prokfinder-refdb`, `--refdb-dir`) are shared by every run. Each database folder has
`COMPLETED-DOWNLOAD.txt` with `format=<n>`; when the format in the folder differs from `PF_REFDB_FORMAT_<DB>`,
search-genomes/download-wf stop and `setup` rebuilds the tables from the raw files already downloaded (nothing is
downloaded again). Each database also writes `NULL-VALUES-<DB>.tsv` (values converted to NULL).

## Conventions (follow them in every change)

- **Language**: code, comments, help and log messages in English; `README.md` in Portuguese.
- **Bash**: `#!/usr/bin/env bash`, `set -euo pipefail`, bash >= 4.4 (the UFV cluster). Variables in `${UPPER}`,
  quoted. Every script sources `lib/common.sh` and uses `log "message|warning|error"` and `die`; temporary files go
  to `make_temp_dir` (removed on exit).
- **Never pipe a large output into `head` or `grep -q`** under `pipefail` (SIGPIPE, exit 141): use `awk 'NR <= n'`
  or bash tests.
- **Command line**: one flag = one criterion; no environment variables for parameters (only CHECKM2DB,
  GTDBTK_DATA_PATH and NCBI_API_KEY, which configure the programs). Help only with `-h`; a command-line error shows
  only the message and which `-h` to read. Filters must be joined with AND/OR.
- **Tables**: TAB-separated, header in UPPER_CASE for Prokfinder columns, empty/invalid values written as `NULL`
  through the single rule `PF_AWK_NULL` / `format_tsv` (never a new ad hoc rule).
- **Reference table change** -> increase `PF_REFDB_FORMAT_<DB>` in `lib/common.sh` and make the module rebuild the
  table from the raw files.
- **Re-execution**: finished work is reused. Per-genome results are kept in cache tables (`pending_for_cache`,
  `update_cache`) keyed by accession and file modification time; a run never recomputes what did not change.
- **Programs from conda environments** are called with `run_in_own_env <program>` (otherwise another environment's
  python3 may be used; this broke BUSCO on the cluster).
- **Resources**: programs use `-p` CPUs and `--memory` GB; never more than a PBS job asked for.
- **Taxonomy is GTDB only** (no NCBI taxonomy filters). NCBI TaxIDs are used only to link BacDive.
- Module headers: line 3 is `# <file> - Prokfinder module <step>: <description>` and a `# Usage:` block
  (tools/build-catalog.py reads them).
- After editing from Windows (`\\wsl.localhost` paths) the executable bit is lost: run
  `chmod +x main-prokfinder.sh modules/*.sh modules/*.py modules/plot-statistics.R docker/*.sh tools/*.py`.

## Testing

There is no automated test suite yet. Before committing:

1. `for f in main-prokfinder.sh lib/common.sh modules/*.sh docker/*.sh; do bash -n "$f"; done` and
   `shellcheck -S warning <changed files>`.
2. Command-line checks need no databases: invalid options, `run-on-cluster --dry-run`, `--check`, `<sub-command> -h`.
3. Behaviour checks use a small fixture: real slices of the reference tables (e.g., the rows of Thermus and
   Pyrococcus of the NCBI, GTDB and BacDive tables) in a folder given with `--refdb-dir`, and fake programs in the
   PATH (shell scripts named `gtdbtk`, `checkm2`, `datasets`... that write the files the module reads).
4. Real runs are made on the UFV cluster (PBS) with `run-on-cluster`; the login node kills commands > 10 min.
5. `python3 tools/build-catalog.py --check` must pass (regenerate the catalog after changing options or modules).

## Git

- `main` is the published branch; development happens on feature branches (currently `restructure-pipeline`) and
  is merged after a real test on the cluster. Ask the maintainer before merging into `main`.
- Pushing to the branch `docker` (or to `main`, or publishing a release) builds and pushes the Docker image
  `jprocha2005/prokfinder`.
- Commit messages: imperative summary line in English, then what changed and why.
- Update `CHANGELOG.md` (`Unreleased`) and `catalog.json` in the same commit as the change (see
  `.claude/skills/update-docs/SKILL.md`).
