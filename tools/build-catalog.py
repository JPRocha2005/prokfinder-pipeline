#!/usr/bin/env python3
"""
build-catalog.py - writes catalog.json, the machine-readable index of Prokfinder (for people and AI agents).

Everything is read from the code itself, so the catalog never has to be edited by hand:
  - version, sub-commands and command-line options (with their section and help text) from main-prokfinder.sh
  - dependencies (DEPS_* lists) and environment variables of the programs
  - modules: step, description and usage from the header comment of each file in modules/
  - shared functions of lib/common.sh (name and the comment line above them)
  - table formats of the reference databases (PF_REFDB_FORMAT_* in lib/common.sh)

Usage (from the root of the repository):
  python3 tools/build-catalog.py           write catalog.json
  python3 tools/build-catalog.py --check   exit 1 if catalog.json is not up to date (nothing is written)
"""

import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAIN = os.path.join(ROOT, "main-prokfinder.sh")
COMMON = os.path.join(ROOT, "lib", "common.sh")
MODULES = os.path.join(ROOT, "modules")
OUTPUT = os.path.join(ROOT, "catalog.json")

# Option lists of main-prokfinder.sh: <variable> -> (section, takes a value)
OPTION_LISTS = {
    "GLOBAL_FLAGS": ("global", False), "GLOBAL_VALUES": ("global", True),
    "SETUP_FLAGS": ("setup", False), "SETUP_VALUES": ("setup", True),
    "SEARCH_FLAGS": ("search-genomes", False), "SEARCH_VALUES": ("search-genomes", True),
    "SEARCH_FILTERS": ("search-genomes", True),
    "DOWNLOAD_FLAGS": ("download-wf", False), "DOWNLOAD_VALUES": ("download-wf", True),
    "CLUSTER_FLAGS": ("run-on-cluster", False), "CLUSTER_VALUES": ("run-on-cluster", True),
}


def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def shell_string(text, name):
    """Value of NAME="..." (first assignment; may span several lines)."""
    match = re.search(r'^%s="(.*?)"' % re.escape(name), text, re.M | re.S)
    return match.group(1) if match else ""


def shell_array(text, name):
    """Items of NAME=( "a" "b" ... )."""
    match = re.search(r"^%s=\((.*?)^\)" % re.escape(name), text, re.M | re.S)
    if not match:
        match = re.search(r"^%s=\((.*?)\)" % re.escape(name), text, re.M)
    return re.findall(r'"([^"]*)"|(\S+)', match.group(1)) if match else []


def options(main):
    """Command-line options: name, aliases, section, whether it takes a value and its help text."""
    # Help lines start with two spaces and the option; a description continued on the next lines is joined to it
    help_lines = []
    for line in main.splitlines():
        if re.match(r"^  -[^ ]", line):   # ("  - text" is a bullet of the help, not an option)
            help_lines.append(line.strip())
        elif help_lines and re.match(r"^\s{20,}\S", line) and not line.strip().startswith(("-", "${", "[[", "#")):
            if re.split(r"\s{2,}", help_lines[-1], maxsplit=1)[1:]:
                help_lines[-1] += " " + line.strip()
        elif line.strip():
            help_lines.append("")
    seen = {}
    for variable, (section, takes_value) in OPTION_LISTS.items():
        for flag in shell_string(main, variable).split():
            seen.setdefault(flag, (section, takes_value, "filter" if variable == "SEARCH_FILTERS" else "option"))
    result = []
    done = set()
    for flag, (section, takes_value, kind) in seen.items():
        if flag in done:
            continue
        description, aliases, argument = "", [flag], ""
        for line in help_lines:
            # "-o, --results-dir <dir>   Output directory ..." -> names and text, separated by 2+ spaces
            names, text = (re.split(r"\s{2,}", line, maxsplit=1) + [""])[:2]
            words = [w.strip(",") for w in names.split() if w.startswith("-")]
            if flag in words:
                aliases = words
                argument = " ".join(w for w in names.split() if w.startswith("<"))
                description = text.strip()
                break
        done.update(aliases)
        result.append({
            "name": max(aliases, key=len), "aliases": sorted(aliases, key=len), "section": section,
            "kind": kind, "takes_value": takes_value, "argument": argument, "description": description,
        })
    return result


def dependencies(main):
    deps = {}
    for variable, group in (("DEPS_SEARCH", "search"), ("DEPS_DOWNLOAD", "download-wf"), ("DEPS_OPTIONAL", "optional")):
        deps[group] = []
        for quoted, bare in shell_array(main, variable):
            program, _, install = (quoted or bare).partition("|")
            deps[group].append({"program": program, "install": install})
    deps["download-wf"].append({"program": "pyrodigal", "install": shell_string(main, "PYRODIGAL_HINT")})
    deps["unix_tools"] = [quoted or bare for quoted, bare in shell_array(main, "DEPS_SYSTEM")]
    deps["environment_variables"] = {
        "CHECKM2DB": "path of the CheckM2 database file (download-wf)",
        "GTDBTK_DATA_PATH": "folder of the GTDB-Tk reference package (ANI and tree classification)",
        "NCBI_API_KEY": "optional NCBI API key (more requests per second)",
    }
    return deps


def modules():
    result = []
    for name in sorted(os.listdir(MODULES)):
        path = os.path.join(MODULES, name)
        if not os.path.isfile(path) or not name.endswith((".sh", ".py", ".R")):
            continue
        lines = read(path).splitlines()
        header = []
        for line in lines[1:]:
            if line.startswith('"""'):
                continue
            if name.endswith(".py"):
                if line.startswith('"""') or line.startswith("import"):
                    break
                header.append(line)
                continue
            if not line.startswith("#"):
                break
            header.append(line[1:].lstrip(" ") if line.startswith("# ") else line[1:])
        title = next((line for line in header if line.strip()), "")
        step, description = "", title.split(" - ", 1)[-1]
        description = re.sub(r"^Prokfinder: ", "", description)
        match = re.match(r"Prokfinder module ([0-9.]+|\([^)]*\))?:?\s*(.*)", description)
        if match:
            step, description = (match.group(1) or "").strip("()"), match.group(2)
        usage = [line.strip() for line in header if re.match(r"\s*(Usage:|[a-z-]+\.(sh|py) )", line)]
        result.append({
            "file": "modules/" + name, "step": step, "description": description.strip(),
            "usage": usage, "language": {"sh": "bash", "py": "python", "R": "R"}[name.rsplit(".", 1)[1]],
        })
    return result


def common_functions(common):
    result = []
    lines = common.splitlines()
    for index, line in enumerate(lines):
        match = re.match(r"^([a-z_][a-z0-9_]*)\(\)\s*\{", line)
        if not match:
            continue
        comment = []
        k = index - 1
        while k >= 0 and lines[k].startswith("#"):
            comment.insert(0, lines[k].lstrip("# ").rstrip())
            k -= 1
        usage = next((c for c in comment if c.startswith("Usage:")), "")
        summary = " ".join(c for c in comment if c and not c.startswith("Usage:"))
        result.append({"name": match.group(1), "summary": summary, "usage": usage})
    return result


def build():
    main = read(MAIN)
    common = read(COMMON)
    return {
        "name": "Prokfinder",
        "version": shell_string(main, "PROKFINDER_VERSION"),
        "description": "Selection, download and curation of prokaryote genomes (NCBI metadata, BacDive phenotypes, "
                       "GTDB taxonomy)",
        "repository": "https://github.com/LBCM-UFV/prokfinder-pipeline",
        "license": "GPL-3.0-or-later",
        "entry_point": "main-prokfinder.sh",
        "generated_by": "tools/build-catalog.py (do not edit by hand)",
        "subcommands": shell_string(main, "SUBCOMMANDS").split(),
        "options": options(main),
        "modules": modules(),
        "common_functions": common_functions(common),
        "reference_database_formats": {
            db: int(value) for db, value in re.findall(r"^PF_REFDB_FORMAT_([A-Z]+)=([0-9]+)", common, re.M)
        },
        "dependencies": dependencies(main),
    }


def main():
    text = json.dumps(build(), indent=2, ensure_ascii=False) + "\n"
    if "--check" in sys.argv[1:]:
        current = read(OUTPUT) if os.path.exists(OUTPUT) else ""
        if current != text:
            print("catalog.json is out of date: run python3 tools/build-catalog.py", file=sys.stderr)
            return 1
        print("catalog.json is up to date")
        return 0
    with open(OUTPUT, "w", encoding="utf-8", newline="\n") as handle:
        handle.write(text)
    print("catalog.json written")
    return 0


if __name__ == "__main__":
    sys.exit(main())
