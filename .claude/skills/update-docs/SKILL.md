---
name: update-docs
description: Update the documentation and indexes of Prokfinder after a code change (catalog.json, CHANGELOG.md, README.md, CLAUDE.md, reference table formats). Use it before every commit that changes options, modules, outputs, dependencies or reference tables, and when preparing a release.
---

# Update the documentation after a code change

Run this checklist before committing a change to `main-prokfinder.sh`, `lib/common.sh`, `modules/` or `docker/`.

## 1. Regenerate the catalog

```bash
python3 tools/build-catalog.py          # rewrites catalog.json from the code
python3 tools/build-catalog.py --check  # must print "catalog.json is up to date"
```

Never edit `catalog.json` by hand. If an option or module is missing from it, fix the source instead: options must
be in the `*_FLAGS`/`*_VALUES`/`SEARCH_FILTERS` lists and in the help of `main-prokfinder.sh` (a line starting with
two spaces and the option); modules need the header line `# <file> - Prokfinder module <step>: <description>`.

## 2. CHANGELOG.md

Add one line per user-visible change under `## [Unreleased]`, in the right group: `Added`, `Changed`, `Removed`,
`Fixed`. Write what the user notices (option, output, behaviour), not the implementation. Internal refactoring
with no visible effect does not go in.

## 3. Reference tables

If a table written by `setup` changed (columns, values, NULL rule): increase `PF_REFDB_FORMAT_<DB>` in
`lib/common.sh`, make sure the module rebuilds the table from the raw files, and add a `Changed` line with the new
format number. The Zenodo packages of that database must be rebuilt.

## 4. README.md (Portuguese)

Update the sections that describe what changed: options tables, steps of `download-wf`, output structure,
dependencies, re-execution. Keep the README in Portuguese and the help/log messages in English.

## 5. CLAUDE.md

Only when the architecture, a convention or the way to test changes (a new module, a new shared rule, a new
folder). Keep it short.

## 6. Before committing

```bash
for f in main-prokfinder.sh lib/common.sh modules/*.sh docker/*.sh; do bash -n "$f"; done
shellcheck -S warning <changed .sh files>
chmod +x main-prokfinder.sh modules/*.sh modules/*.py modules/plot-statistics.R docker/*.sh tools/*.py
```

Commit the code, `catalog.json`, `CHANGELOG.md` and the documentation together.

## At a release

1. Rename `## [Unreleased]` to `## [X.Y.Z] - YYYY-MM-DD` and open a new empty `## [Unreleased]`.
2. Set `PROKFINDER_VERSION` in `main-prokfinder.sh` and `version`/`date-released` in `CITATION.cff`.
3. Regenerate `catalog.json`, merge into `main` and publish the GitHub release `vX.Y.Z` (it builds the Docker image).
