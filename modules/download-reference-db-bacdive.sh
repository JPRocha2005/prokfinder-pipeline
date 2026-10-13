#!/usr/bin/env bash
#
# download-reference-db-bacdive.sh - Prokfinder module 0.2: BacDive reference database
#
# BactoTraits Script #1 (SCRIPTS/1__bacdive_fetch().R) rewritten in Bash, using the BacDive API exclusively.
# The jq programs are embedded below, so no extra files are needed. Every table written is TAB-separated.
#
#   1. downloads the list of BacDive IDs (advanced search CSV -> TSV with header)
#   2. downloads the JSON records of every strain through the API, in batches (interrupted runs are resumed)
#   3. extracts the NCBI TaxIDs and 29 traits of each strain into one TSV (header without spaces, empty -> NULL)
#
# Usage:
#    download-reference-db-bacdive.sh -o <bacdive_dir>
#    download-reference-db-bacdive.sh -o <bacdive_dir> --limit 5000
#    download-reference-db-bacdive.sh -o <bacdive_dir> --skip-download   # only re-extract the traits
#
# Output inside <bacdive_dir>:
#    BACDIVE-DOWNLOAD-BRUTE/advsearch_bacdive.tsv       (list of BacDive IDs, TSV with header)
#    BACDIVE-DOWNLOAD-BRUTE/batch_*.jsonl               (raw API results in JSON, one strain per line)
#    BACDIVE-DOWNLOAD-BRUTE/list_download_errors.tsv    (IDs that failed after all retry attempts)
#    BACDIVE-FORMATTED-TSV/BACDIVE_DATA_FORMATTED.tsv   (one row per strain: bacdive_id, ncbi_tax_id and one
#                                                       column per trait)
#    BACDIVE-FORMATTED-TSV/BACDIVE_DATA_UNFORMATTED.tsv (OPTIONAL, --unformatted: long dump
#                                                       bacdive_id, L2, L3, L4, L5, value)
#    COMPLETED-DOWNLOAD.txt                             (checkpoint: written only when every batch was downloaded)
#
# Format of a trait cell, e.g. temperature_preference = "20:35:40;37":
#    ";" separates references (each reference is a different publication or culture collection record)
#    ":" separates values of the same reference (a range "20-40" is kept as its two endpoints, 20:40)
#    ncbi_tax_id holds every distinct NCBI TaxID of the strain (e.g. species and strain level), joined with ";"
#  *_opt        values recorded as "optimum"
#  *_preference every growth-positive value, whatever its label (growth, optimum, minimum or maximum)
#  *_range      only values recorded as minimum or maximum (so "10:40" means minimum 10 and maximum 40)
#  *_delta      max - min of all growth-positive values (empty when there are fewer than two distinct values)

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../lib/common.sh
source "${SCRIPT_DIR}/../lib/common.sh"

API_URL="${BACDIVE_API_URL:-https://api.bacdive.dsmz.de/v2}"
BACDIVE_CSV_URL="https://bacdive.dsmz.de/advsearch/csv"

BATCH_SIZE=100           # Request up to 100 strains per single API call
MAX_RESTARTS=3           # Retry up to 3 times per batch
RETRY_WAIT=5             # Wait 5 seconds before each retry

OUT_DIR=""
IDS_FILE=""
LIMIT=""
SKIP_DOWNLOAD=0
MAKE_UNFORMATTED=0
FORCE=0

# The 29 BacDive-derived traits, in output column order
TRAITS=(
    gram_stain cell_length cell_width cell_shape motility flagellum_arrangement
    pigmentation_production colony_color forms_multicellular_complex
    temperature_opt temperature_preference temperature_range temperature_delta
    pH_opt pH_preference pH_range pH_delta
    spore_formation
    halophily_opt halophily_preference halophily_range halophily_delta
    oxygen_tolerance antibiotic_resistance_frequency
    energy_source electron_donor carbon_source
    biosafety_level gc_content
)

usage() {
    cat <<EOF
Usage: $(basename "$0") -o|--out-dir DIR [--ids-file FILE] [--limit N] [--retry-wait S] [--batch-size N]
          [--skip-download] [--unformatted] [--force]

Required:
  -o, --out-dir DIR    Output directory (Prokfinder/INPUT/REFERENCE-METADATA/BACDIVE).

Options:
  --ids-file FILE      Alternative TSV/CSV file whose first column holds BacDive IDs.
                       Default: downloaded to BACDIVE-DOWNLOAD-BRUTE/advsearch_bacdive.tsv
  --limit N            Process only the first N IDs (useful for testing)
  --retry-wait S       Wait S seconds between retries (default: 5)
  --batch-size N       Strains requested per API request (default: 100)
  --skip-download      Skip the API download; only rebuild the tables from the batch files already downloaded
  --unformatted        Also write BACDIVE_DATA_UNFORMATTED.tsv (generic long dump of every field; large)
  --force              Remove the previous download of this module and start again
EOF
    exit "${1:-1}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        -o|--out-dir) OUT_DIR="${2%/}"; shift 2 ;;
        --ids-file|--ids-csv) IDS_FILE="$2"; shift 2 ;;
        --limit) LIMIT="$2"; shift 2 ;;
        --retry-wait) RETRY_WAIT="$2"; shift 2 ;;
        --batch-size) BATCH_SIZE="$2"; shift 2 ;;
        --skip-download) SKIP_DOWNLOAD=1; shift ;;
        --unformatted) MAKE_UNFORMATTED=1; shift ;;
        --force) FORCE=1; shift ;;
        -h|--help) usage 0 ;;
        *) log "error" "Unknown option: $1"; usage ;;
    esac
done

if [ -z "$OUT_DIR" ]; then
    log "error" "Mandatory option -o|--out-dir is missing."
    usage
fi

DOWNLOAD_DIR="$(pf_bacdive_brute_dir "$OUT_DIR")"
FORMATTED="$(pf_bacdive_table "$OUT_DIR")"
FORMATTED_DIR="$(dirname "$FORMATTED")"
COMPLETED_MARKER="$(pf_completed_marker "$OUT_DIR")"
IDS_FILE="${IDS_FILE:-$DOWNLOAD_DIR/advsearch_bacdive.tsv}"
ERROR_FILE="$DOWNLOAD_DIR/list_download_errors.tsv"

mkdir -p "$OUT_DIR"
setup_module_log "$OUT_DIR/download-reference-db-bacdive.log"
require_programs curl jq gawk

if [ "$FORCE" -eq 1 ]; then
    log "message" "Option --force: removing the previous BacDive download in $OUT_DIR"
    rm -rf "$DOWNLOAD_DIR" "$FORMATTED_DIR"
    rm -f "$COMPLETED_MARKER"
fi

mkdir -p "$DOWNLOAD_DIR" "$FORMATTED_DIR"

# ===========================================================================
# EMBEDDED jq PROGRAM 1 — extraction of the NCBI tax id and the 29 traits (one TSV row per strain)
#
# Input : one JSON object per line, {"id": "<bacdive_id>", "data": {...}}
# Arg   : $cols = JSON array with the trait names, in output order
# Output: bacdive_id <TAB> ncbi_tax_id <TAB> trait1 <TAB> trait2 ...
#
# Cell format (same idea as the original R script):
#   - values that come from the SAME reference (@ref)   are joined with ":"
#   - values that come from DIFFERENT references        are joined with ";"
#   - everything is alphabetically/numerically ordered ("alpha_order")
#   - ncbi_tax_id is not a trait: it holds every distinct NCBI taxonomy ID of the
#     strain (e.g. species and strain level), joined with ";"
#
# If BacDive renames a field, only the key names in the trait definitions at
# the bottom of this program need to be edited.
# ===========================================================================
JQ_TRAITS=$(cat <<'JQ_EOF'
# ---------- generic helpers -------------------------------------------------
def arr: if . == null then [] elif type == "array" then . else [.] end;

# first non-null, non-empty value among the alternatives
def pick(fs): [fs] | map(select(. != null and . != "")) | first // null;

def ref: (.["@ref"] // "NA") | tostring;

# records of a subsection (BacDive gives either one object or a list of objects)
def recs(sec; sub):
    ((.[sec] // {}) | if type == "object" then .[sub] else null end)
    | arr | map(select(type == "object"));

# numbers inside a string: "37" -> [37]; "5-30" -> [5,30]; "1.5-3 µm" -> [1.5,3]
def nums:
    if (type == "string" or type == "number") then tostring else "" end
    | gsub("[^0-9.\\-]"; "")
    | if test("^-?[0-9]+(\\.[0-9]+)?-[0-9]+(\\.[0-9]+)?$")
      then capture("^(?<a>-?[0-9]+(\\.[0-9]+)?)-(?<b>[0-9]+(\\.[0-9]+)?)$") | [.a, .b | tonumber]
      elif test("^-?[0-9]+(\\.[0-9]+)?$") then [tonumber]
      else [] end;

# list of {ref, v}  ->  "v1:v2;v3"
def cell:
    if length == 0 then ""
    else group_by(.ref)
         | map(map(.v) | unique | map(tostring) | join(":"))
         | map(select(. != "")) | sort | join(";")
    end;

def nominal(sec; sub; vf; mapper):
    [ recs(sec; sub)[] | ref as $r
      | (vf // empty | tostring | ascii_downcase | mapper) as $v
      | select($v != "")
      | {ref: $r, v: $v} ];

def numeric(sec; sub; vf):
    [ recs(sec; sub)[] | ref as $r | (vf // empty | nums[]) as $v | {ref: $r, v: $v} ];

# ---------- NCBI taxonomy ID (General section) ------------------------------
# BacDive gives either one object or a list of objects, e.g.
#   "NCBI tax id": {"NCBI tax id": 2926319, "Matching level": "species"}
# All distinct IDs are kept (a strain can be linked at species and at strain level),
# sorted numerically and joined with ";". Empty when the strain has no NCBI tax id.
def ncbi_taxids:
    [ ((.General // {}) | if type == "object" then .["NCBI tax id"] else null end)
      | arr[]
      | (if type == "object" then .["NCBI tax id"] else . end)
      | select(. != null and . != "")
      | tostring ]
    | unique
    | sort_by(tonumber? // 0)
    | join(";");

# ---------- value harmonisation (modalities described in the BactoTraits paper)
def yesno:
    if test("^(yes|true|positive|\\+)$") then "yes"
    elif test("^(no|false|negative|-)$") then "no"
    else empty end;

def shape:
    if   test("oval|ellips|ovoid|egg")                       then "ovoid"
    elif test("spiral|helical|spirill")                      then "spiral"
    elif test("cocc|sphere|spheric")                         then "coccus"
    elif test("vibrio|crescent|curved|comma|bent")           then "curved"
    elif test("rod|bacill")                                  then "rod"
    elif test("pleomorph")                                   then "pleomorphic"
    elif test("filament|hypha|thread")                       then "filament"
    else "other" end;

def flagellum:
    if   test("peritrich")                                                    then "peritrichous"
    elif test("gliding")                                                      then "gliding"
    elif test("polar|monotrich|lophotrich|amphitrich|tuft|polytrich")         then "polar"
    else empty end;

def colour:
    . as $s
    | ( (if ($s | test("colou?rless|translucent|transparent|hyaline")) then "colorless" else empty end),
        (if ($s | test("\\b(cream|beige|ivory|buff)"))                 then "cream"     else empty end),
        (if ($s | test("\\bgr[ae]y"))                                  then "grey"      else empty end),
        ( ("yellow","white","brown","pink","orange","red","green","blue") as $c
          | if ($s | test("\\b" + $c)) then $c else empty end ) );

def oxygen:
    if   test("microaero")   then "microaerophile"
    elif test("facultative") then "facultative"
    elif test("anaerob")     then "anaerobe"
    elif test("aerob")       then "aerobe"
    else empty end;

# nutrition type, e.g. "chemoorganoheterotroph" -> energy / electron donor / carbon source
def energy:  (if test("chemo")   then "chemotroph"   else empty end), (if test("photo")  then "phototroph"   else empty end);
def donor:   (if test("organo")  then "organotroph"  else empty end), (if test("litho")  then "lithotroph"   else empty end);
def carbon:  (if test("hetero")  then "heterotroph"  else empty end), (if test("autotro") then "autotroph"   else empty end);

# all string values of a record except @ref
def strvals: to_entries | map(select(.key != "@ref" and (.value | type) == "string")) | map(.value)[];

# ---------- temperature / pH / salinity ------------------------------------
# keeps only growth-positive records that carry at least one number
def envrecs(sec; sub; valf):
    [ recs(sec; sub)[]
      | ref as $r
      | ((pick(.growth, .ability) // "") | ascii_downcase) as $ab
      | ((pick(.type, .["tested relation"]) // "") | ascii_downcase) as $ty
      | select($ab != "negative" and $ab != "no" and ($ty | test("no growth") | not))
      | {ref: $r, type: $ty, nums: ((valf // "") | nums)}
      | select(.nums | length > 0) ];

def env_opt(rs):  [ rs[] | select(.type == "optimum") | {ref, v: .nums[]} ];
def env_rng(rs):  [ rs[] | select(.type == "minimum" or .type == "maximum") | {ref, v: .nums[]} ];
def env_all(rs):  [ rs[] | {ref, v: .nums[]} ];
# max - min of all growth-positive values; left empty when the strain has fewer than
# two distinct values (a single test point says nothing about tolerance).
# To get 0 instead, change "< 2" to "== 0".
def env_delta(rs):
    ([rs[].nums[]] | unique) as $a
    | if ($a | length) < 2 then ""
      else (($a | last) - ($a | first)) | ((. * 1000 | round) / 1000 | tostring) end;

# NaCl concentration in %, skipping molar / g per litre values
def salt_pct:
    select(((.salt // "NaCl") | tostring | ascii_downcase) == "nacl")
    | (.concentration // .["salt concentration"] // .["Salt conc."])
    | select(. != null) | tostring
    | select(test("/\\s*[lL]\\b|mol|\\bm?M\\b") | not);

# ratio of resistances found / molecules tested, computed per reference
def antibiotic_ratio:
    [ recs("Physiology and metabolism"; "antibiotic resistance")[]
      | {ref: ref,
         res:  ((.["is resistant"]    // "") | tostring | ascii_downcase),
         sens: ((.["is sensitive"]    // "") | tostring | ascii_downcase),
         inter:((.["is intermediate"] // "") | tostring | ascii_downcase)}
      | select(.res != "" or .sens != "" or .inter != "") ]
    | group_by(.ref)
    | map(((map(select(.res == "yes")) | length) / length) | ((. * 1000 | round) / 1000))
    | sort | map(tostring) | join(";");

# ---------- one row per strain ---------------------------------------------
.id as $id
# The NCBI tax id is computed outside the traits "try": a problem in a trait must not lose the ID that links the strain to NCBI
| (try ((.data // {}) | ncbi_taxids) catch "") as $tax
| (try (
    (.data // {}) as $d
    | ($d | envrecs("Culture and growth conditions"; "culture temp"; .temperature)) as $T
    | ($d | envrecs("Culture and growth conditions"; "culture pH";   .pH))          as $P
    | ($d | envrecs("Physiology and metabolism";     "halophily";    salt_pct))    as $H
    | ($d | {
        gram_stain:  (nominal("Morphology"; "cell morphology"; .["gram stain"]; select(. == "positive" or . == "negative")) | cell),
        cell_length: (numeric("Morphology"; "cell morphology"; pick(.["cell length"], .["cell len"])) | cell),
        cell_width:  (numeric("Morphology"; "cell morphology"; pick(.["cell width"])) | cell),
        cell_shape:  (nominal("Morphology"; "cell morphology"; .["cell shape"]; shape) | cell),
        motility:    (nominal("Morphology"; "cell morphology"; .motility; yesno) | cell),
        flagellum_arrangement: (nominal("Morphology"; "cell morphology"; .["flagellum arrangement"]; flagellum) | cell),
        pigmentation_production: (nominal("Morphology"; "pigmentation"; pick(.production, .["pigment production"]); yesno) | cell),
        colony_color: (nominal("Morphology"; "colony morphology"; pick(.["colony color"], .color, .["colony colour"]); colour) | cell),
        forms_multicellular_complex: (nominal("Morphology"; "multicellular morphology";
                                        pick(.["forms multicellular complex"], .ability, .["multicellular complex forming ability"]); yesno) | cell),
        spore_formation: (nominal("Physiology and metabolism"; "spore formation"; pick(.["spore formation"], .ability); yesno) | cell),
        oxygen_tolerance: (nominal("Physiology and metabolism"; "oxygen tolerance"; pick(.["oxygen tolerance"], .["oxygen tol"]); oxygen) | cell),
        energy_source:    (nominal("Physiology and metabolism"; "nutrition type"; strvals; energy) | cell),
        electron_donor:  (nominal("Physiology and metabolism"; "nutrition type"; strvals; donor)  | cell),
        carbon_source:    (nominal("Physiology and metabolism"; "nutrition type"; strvals; carbon) | cell),
        biosafety_level:  (nominal("Interaction and safety"; "risk assessment"; .["biosafety level"]; select(test("^[0-9]+$"))) | cell),
        gc_content:       (numeric("Sequence information"; "GC content"; pick(.["GC-content"], .["GC content"])) | cell),
        antibiotic_resistance_frequency: antibiotic_ratio
      }
      + {
        temperature_opt:        (env_opt($T)  | cell),
        temperature_preference: (env_all($T)  | cell),
        temperature_range:      (env_rng($T)  | cell),
        temperature_delta:      env_delta($T),
        pH_opt:                 (env_opt($P)  | cell),
        pH_preference:          (env_all($P)  | cell),
        pH_range:               (env_rng($P)  | cell),
        pH_delta:               env_delta($P),
        halophily_opt:          (env_opt($H)  | cell),
        halophily_preference:   (env_all($H)  | cell),
        halophily_range:        (env_rng($H)  | cell),
        halophily_delta:        env_delta($H)
      })
    ) catch ("strain \($id): \(.)" | debug | {})) as $t
| [$id, $tax] + [ $cols[] as $c | ($t[$c] // "") ]
| @tsv
JQ_EOF
)

# ===========================================================================
# EMBEDDED jq PROGRAM 2 — generic long dump (only used with --unformatted)
#   bacdive_id, L2 (section), L3, L4, L5 (deeper levels joined with "|"), value
# ===========================================================================
JQ_MELT=$(cat <<'JQ_EOF'
.id as $id
| (.data // empty)
| tostream
| select(length == 2)
| select((.[1] | type) as $t | $t != "array" and $t != "object" and $t != "null")
| (.[0] | map(select(type == "string"))) as $p
| .[1] as $v
| select($v != "")
| [ $id, ($p[0] // ""), ($p[1] // ""), ($p[2] // ""), ($p[3:] | join("|")), $v ]
| @tsv
JQ_EOF
)

# ---------------------------------------------------------------------------
# 0/3 — DOWNLOAD AND CLEANUP OF ADVANCED SEARCH CSV  ->  TSV (with header)
# ---------------------------------------------------------------------------
download_advsearch_tsv() {
    local target="$1"
    log "message" "Downloading advanced search CSV from BacDive: $BACDIVE_CSV_URL"

    local tmp_raw
    tmp_raw="$(mktemp)"

    if ! curl -sSL --retry 3 "$BACDIVE_CSV_URL" -o "$tmp_raw"; then
        rm -f "$tmp_raw"
        die "Could not download the CSV file from BacDive."
    fi

    # 1. tail -n +4 : skip the first 3 notice lines from BacDive
    # 2. gawk FPAT  : parse the (quoted) CSV correctly and re-emit it as TSV. The first row whose first field
    #                 is not a number is kept as the header (spaces -> "_"; "bacdive_id" if it has no name);
    #                 after it, only rows whose first field is a numeric BacDive ID are kept
    tail -n +4 "$tmp_raw" | tr -d '\r' | gawk -v FPAT='([^,]*)|("([^"]|"")*")' -v OFS='\t' '
        function clean(f) {
            if (f ~ /^".*"$/) { f = substr(f, 2, length(f) - 2); gsub(/""/, "\"", f) }
            gsub(/[\t\n]/, " ", f)
            return f
        }
        {
            id = clean($1)
            if (id !~ /^[0-9]+$/) {
                if (header_done || NF < 2) next
                line = ""
                for (i = 1; i <= NF; i++) {
                    h = clean($i)
                    gsub(/^[[:space:]]+|[[:space:]]+$/, "", h)
                    gsub(/[[:space:]]+/, "_", h)
                    if (i == 1 && (h == "" || tolower(h) == "id")) h = "bacdive_id"
                    line = line (i > 1 ? OFS : "") h
                }
                print line
                header_done = 1
                next
            }
            if (!header_done) { print "bacdive_id"; header_done = 1 }
            line = ""
            for (i = 1; i <= NF; i++) line = line (i > 1 ? OFS : "") clean($i)
            print line
        }
    ' > "$target"

    rm -f "$tmp_raw"
    log "message" "Clean TSV successfully saved to: $target ($(count_records "$target") IDs)"
}

# ---------------------------------------------------------------------------
# 1/3 — ID LIST EXTRACTION (first column of a TSV or CSV; non-numeric lines, such as the header, are skipped)
# ---------------------------------------------------------------------------
get_id_list() {
    gawk -F'[\t,]' '{ gsub(/"/, "", $1); gsub(/\r/, "", $1); if ($1 ~ /^[0-9]+$/) print $1 }' "$1"
}

# ---------------------------------------------------------------------------
# 2/3 — DATA RETRIEVAL (Multi-ID endpoint: /fetch/id1;id2;id3...)
# ---------------------------------------------------------------------------
fetch_batch() {
    local joined_ids="$1"
    local attempt=0
    local body http_code

    while [ "$attempt" -le "$MAX_RESTARTS" ]; do
        if [ "$attempt" -gt 0 ]; then
            log "message" "Batch transfer error, waiting ${RETRY_WAIT}s before retry $attempt/$MAX_RESTARTS"
            sleep "$RETRY_WAIT"
        fi

        body="$(curl -s -w $'\n%{http_code}' --max-time 120 "$API_URL/fetch/$joined_ids" 2>/dev/null)"
        http_code="${body##*$'\n'}"
        body="${body%$'\n'*}"

        if [ "$http_code" = "200" ] && echo "$body" | jq -e . >/dev/null 2>&1; then
            printf '%s' "$body"
            return 0
        fi

        attempt=$((attempt + 1))
    done

    return 1
}

# Downloads every batch that is not on disk yet. Returns 1 if some batch failed (the next run retries it).
download_bacdive() {
    # Download and process the advanced search CSV if no ID file is present locally
    if [ ! -s "$IDS_FILE" ]; then
        download_advsearch_tsv "$IDS_FILE"
    fi

    log "message" "Reading ID list from $IDS_FILE"

    local all_ids=()
    mapfile -t all_ids < <(get_id_list "$IDS_FILE" | { if [ -n "$LIMIT" ]; then head -n "$LIMIT"; else cat; fi; })

    local total="${#all_ids[@]}"
    log "message" "Total IDs to download: $total"
    if [ "$total" -eq 0 ]; then
        die "No BacDive ID found in $IDS_FILE"
    fi

    # Batch files already on disk are kept on purpose (resume). Only the failures of this run are listed.
    echo "bacdive_id" > "$ERROR_FILE"

    local processed=0 batch_n=0 ok=0 failed_batches=0
    local batch_ids joined_ids json current_count batch_file

    while [ "$processed" -lt "$total" ]; do
        batch_ids=("${all_ids[@]:processed:BATCH_SIZE}")

        if [ "${#batch_ids[@]}" -eq 0 ]; then
            break
        fi

        batch_file="$DOWNLOAD_DIR/batch_$(printf '%05d' "$batch_n").jsonl"

        # Resume logic: Check if the batch file already exists and is non-empty
        if [ -s "$batch_file" ]; then
            current_count=$(wc -l < "$batch_file")
            ok=$((ok + current_count))
        else
            # File is missing or empty; perform API request
            : > "$batch_file"

            # Join IDs with semicolons
            joined_ids="$(printf '%s;' "${batch_ids[@]}")"
            joined_ids="${joined_ids%;}"

            if json="$(fetch_batch "$joined_ids")"; then
                # Parse results and append each strain as a JSON line (.results[id])
                jq -c '.results | to_entries[] | {id: .key, data: .value}' <<< "$json" >> "$batch_file"
                current_count=$(wc -l < "$batch_file")
                ok=$((ok + current_count))
            else
                log "warning" "Batch starting at index $((processed + 1)) failed after $MAX_RESTARTS retries — logging IDs to $ERROR_FILE"
                printf '%s\n' "${batch_ids[@]}" >> "$ERROR_FILE"
                failed_batches=$((failed_batches + 1))
            fi

            # Print progress only when a batch was actually fetched
            log "message" "Progress: $((processed + ${#batch_ids[@]}))/$total requested ($ok strains retrieved)"
        fi

        processed=$((processed + ${#batch_ids[@]}))
        batch_n=$((batch_n + 1))
    done

    log "message" "Download finished: $ok/$total strains retrieved successfully."

    if [ "$failed_batches" -gt 0 ]; then
        log "warning" "$failed_batches batch(es) failed. Run the module again to download only the missing batches."
        return 1
    fi
    return 0
}

# ---------------------------------------------------------------------------
# 3/3 — TRAIT EXTRACTION -> BACDIVE_DATA_FORMATTED.tsv
# ---------------------------------------------------------------------------
list_batch_files() {
    shopt -s nullglob
    BATCH_FILES=("$DOWNLOAD_DIR"/*.jsonl)
    shopt -u nullglob
    if [ "${#BATCH_FILES[@]}" -eq 0 ]; then
        log "error" "No batch files (*.jsonl) found in $DOWNLOAD_DIR"
        return 1
    fi
}

melt_batches_to_unformatted() {
    log "message" "Writing generic long dump (bacdive_id, L2, L3, L4, L5, value)..."
    local out="$FORMATTED_DIR/BACDIVE_DATA_UNFORMATTED.tsv"
    printf 'bacdive_id\tL2\tL3\tL4\tL5\tvalue\n' > "$out"

    local f
    for f in "${BATCH_FILES[@]}"; do
        jq -r "$JQ_MELT" "$f" >> "$out" || log "warning" "jq failed on $f"
    done
    log "message" "Long format saved to $out ($(count_records "$out") rows)."
}

format_bacdive() {
    log "message" "Extracting the NCBI tax id and ${#TRAITS[@]} traits from the BacDive JSON batches..."

    local cols_json tmp
    cols_json="$(printf '%s\n' "${TRAITS[@]}" | jq -R . | jq -sc .)"
    tmp="${FORMATTED}.tmp"

    # header: identifiers first (bacdive_id, ncbi_tax_id), then the traits
    { printf 'bacdive_id\tncbi_tax_id'; printf '\t%s' "${TRAITS[@]}"; printf '\n'; } > "$tmp"

    local f
    for f in "${BATCH_FILES[@]}"; do
        # stderr is NOT discarded: per-strain problems are reported by jq's debug()
        if ! jq -r --argjson cols "$cols_json" "$JQ_TRAITS" "$f" >> "$tmp"; then
            rm -f "$tmp"
            die "jq failed on $f"
        fi
    done

    local nstrains
    nstrains=$(count_records "$tmp")
    if [ "$nstrains" -le 0 ]; then
        rm -f "$tmp"
        log "error" "No strains were written to $FORMATTED"
        return 1
    fi

    # Fill summary: how many strains have data for each column (spots wrong key names quickly)
    log "message" "Column fill summary (strains with data / total):"
    gawk -F'\t' '
        NR == 1 { n = NF; for (i = 1; i <= NF; i++) h[i] = $i; next }
        { rows++; for (i = 2; i <= n; i++) if ($i != "") c[i]++ }
        END {
            for (i = 2; i <= n; i++) {
                printf "    %-34s %d/%d%s\n", h[i], c[i] + 0, rows, (c[i] + 0 == 0 ? "   <-- EMPTY: check key names" : "")
            }
        }
    ' "$tmp" >&2

    # Header without spaces and blank cells -> NULL. Done AFTER the fill summary above (which counts empty cells).
    format_tsv 0 < "$tmp" > "$FORMATTED"
    rm -f "$tmp"
    log "message" "Traits saved to $FORMATTED ($nstrains strains; blank cells replaced with NULL)."
}

# ---------------------------------------------------------------------------
main() {
    local download_ok=1

    if [ "$SKIP_DOWNLOAD" -eq 0 ]; then
        download_bacdive || download_ok=0
    else
        log "message" "Skipping download (--skip-download): using existing batches in $DOWNLOAD_DIR"
    fi

    list_batch_files || exit 1
    if [ "$MAKE_UNFORMATTED" -eq 1 ]; then
        melt_batches_to_unformatted || exit 1
    fi
    format_bacdive || exit 1

    # The checkpoint is only written when every batch was downloaded: otherwise Prokfinder calls this module
    # again in its next run, which downloads only the missing batches
    if [ "$download_ok" -eq 1 ]; then
        {
            echo "BacDive reference data downloaded on $(date +"%Y-%m-%d %H:%M:%S")"
            echo "api=$API_URL"
            echo "strains=$(count_records "$FORMATTED")"
        } > "$COMPLETED_MARKER"
        log "message" "BacDive reference data download completed successfully!"
    else
        log "warning" "BacDive table written with the strains downloaded so far; some batches are still missing."
    fi
}

main
