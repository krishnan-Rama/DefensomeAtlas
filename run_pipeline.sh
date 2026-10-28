#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# defensome: end-to-end pipeline for one directory of proteomes.
#
#   bash run_pipeline.sh --proteomes DIR --pfam Pfam-A.hmm --out DIR [options]
#
# Runs, in order: inspect -> collapse -> scan -> annotate -> qc -> extract ->
# domains -> cyp -> trees -> compare -> report -> dashboard.
#
# The first run stops after `inspect` so you can check the sample sheet. That
# is deliberate: every silent failure this pipeline has ever had came from an
# ID or grouping assumption nobody looked at.
# ---------------------------------------------------------------------------
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# Find defensome.py whether this script sits beside it or one level above it
# (unpacking the tarball inside an existing checkout leaves it nested).
DEF=""
for c in "$HERE/defensome.py" "$HERE/defensome/defensome.py" "$HERE/../defensome.py"; do
  [[ -f "$c" ]] && { DEF=$(cd "$(dirname "$c")" && pwd)/defensome.py; break; }
done
if [[ -z "$DEF" ]]; then
  echo "ERROR: defensome.py not found next to $0 or in ./defensome/." >&2
  echo "  The tarball unpacks to a defensome/ directory; run the copy inside it:" >&2
  echo "      bash defensome/run_pipeline.sh ..." >&2
  echo "  or move its contents up one level." >&2
  exit 1
fi

PROT=""; PFAM=""; OUT=""; SHEET=""; META=""; GROUP=""; CYPREFS=""; TREE=""
THREADS=8; REF="DToL"; STOP=""; DROP_UTR=""; SKIP_INSPECT=0; EXCLUDE=""
DOWNLOAD=0; FULL_ARCH=0

usage() { sed -n '2,14p' "$0"; cat <<'U'

Required:
  --out DIR           output directory
  --proteomes DIR     first run only (inspect); ignored once --samplesheet is given
  --pfam FILE         Pfam-A.hmm or Pfam-A.hmm.gz; without it the prebuilt
                      defensome HMMs are downloaded (a few MB)

Optional:
  --samplesheet FILE  use this sheet instead of generating one (skips inspect)
  --metadata FILE     TSV with a trait column for grouping
  --group-by COL      trait column name
  --reference METHOD  method to measure recovery against (default DToL)
  --cyp-refs FILE     clan-labelled CYP reference FASTA
  --species-tree FILE Newick species tree to embed in the dashboard
  --drop-utrorf       discard EvidentialGene utrorf entries when collapsing
  --exclude LIST      comma-separated species or sample_ids to leave out of compare
  --download-pfam     fetch the current Pfam-A from EBI into <out>/db (about 300 MB)
  --full-architecture also report every non-defensome domain on defensome proteins
                      (needs the full Pfam-A; adds minutes, not hours)
  --threads N         default 8
  --stop-after STEP   stop after inspect|collapse|scan|annotate|compare
U
exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --proteomes) PROT=$2; shift 2;;
    --pfam) PFAM=$2; shift 2;;
    --out) OUT=$2; shift 2;;
    --samplesheet) SHEET=$2; SKIP_INSPECT=1; shift 2;;
    --metadata) META=$2; shift 2;;
    --group-by) GROUP=$2; shift 2;;
    --reference) REF=$2; shift 2;;
    --cyp-refs) CYPREFS=$2; shift 2;;
    --species-tree) TREE=$2; shift 2;;
    --drop-utrorf) DROP_UTR="--drop-utrorf"; shift;;
    --exclude) EXCLUDE=$2; shift 2;;
    --download-pfam) DOWNLOAD=1; shift;;
    --full-architecture) FULL_ARCH=1; shift;;
    --threads) THREADS=$2; shift 2;;
    --stop-after) STOP=$2; shift 2;;
    -h|--help) usage;;
    *) echo "unknown option: $1" >&2; usage;;
  esac
done
[[ -n "$OUT" ]] || usage
if [[ "$SKIP_INSPECT" -eq 0 ]]; then
  [[ -n "$PROT" ]] || { echo "ERROR: --proteomes is required for the first (inspect) run" >&2; exit 1; }
  [[ -d "$PROT" ]] || { echo "ERROR: not a directory: $PROT" >&2; exit 1; }
  PROT=$(cd "$PROT" && pwd)
else
  [[ -f "$SHEET" ]] || { echo "ERROR: sample sheet not found: $SHEET" >&2; exit 1; }
  # --proteomes is ignored here: the sheet carries absolute FASTA paths
fi
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd)

# Resolve a file argument to an absolute path, failing if it is not there.
# `readlink -f` returns an EMPTY STRING when a parent directory is missing, so
# canonicalising before checking silently blanked the argument: a placeholder
# --pfam /path/to/Pfam-A.hmm.gz then disabled the scan, skipped the HMMER
# preflight and surfaced four steps later as "no domtblout files".
abs_file() {   # abs_file <label> <path>
  local label=$1 p=$2
  [[ -n "$p" ]] || return 0
  if [[ ! -f "$p" ]]; then
    echo "ERROR: $label not found: $p" >&2
    case "$p" in /path/to/*|*example*|*'<'*) echo "       That looks like a placeholder; give the real path." >&2;; esac
    exit 1
  fi
  ( cd "$(dirname "$p")" && printf '%s/%s\n' "$(pwd)" "$(basename "$p")" )
}
PFAM=$(abs_file "--pfam" "$PFAM") || exit 1
META=$(abs_file "--metadata" "$META") || exit 1
CYPREFS=$(abs_file "--cyp-refs" "$CYPREFS") || exit 1
TREE=$(abs_file "--species-tree" "$TREE") || exit 1
: "${SHEET:=$OUT/samplesheet.tsv}"

FAILED=()
# step: an independent leaf. A failure is recorded and the run continues.
step() {
  echo; echo "########## $1  $(date +%H:%M:%S) ##########"
  if "${@:2}"; then echo "---------- $1 OK"; else
    echo "---------- $1 FAILED (rc=$?)"; FAILED+=("$1"); fi
}
# must: a link in the core chain. Everything downstream depends on it, so a
# failure stops the run rather than producing a cascade of follow-on errors.
must() {
  echo; echo "########## $1  $(date +%H:%M:%S) ##########"
  if "${@:2}"; then echo "---------- $1 OK"; else
    rc=$?
    echo "---------- $1 FAILED (rc=$rc)"
    echo; echo "STOPPING: '$1' is required by every later step."
    echo "  python3 $DEF doctor    shows what is missing"
    exit "$rc"
  fi
}
done_after() { [[ -n "$STOP" && "$STOP" == "$1" ]]; }

# 1. inspect ---------------------------------------------------------------
if [[ "$SKIP_INSPECT" -eq 0 ]]; then
  must inspect python3 "$DEF" inspect --proteomes "$PROT" --out "$OUT" --samplesheet "$SHEET"
  cat <<MSG

========================================================================
A draft sample sheet is at:
  $SHEET

CHECK IT before continuing. The species and method columns are guessed
from file names; gene_regex is guessed from the headers. Everything
downstream trusts this file.

When it is right, rerun with:
  bash run_pipeline.sh --proteomes $PROT --pfam <Pfam-A.hmm> \\
       --out $OUT --samplesheet $SHEET
========================================================================
MSG
  exit 0
fi
done_after inspect && exit 0

# preflight: find out now, not thirty samples in, that a tool is missing
REQ="core,scan"
[[ -n "$CYPREFS" ]] && REQ="$REQ,cyp"
echo "########## preflight ##########"
python3 "$DEF" doctor --require "$REQ" || {
  echo; echo "STOPPING before any work: fix the missing tools above, then rerun."
  exit 1; }

# 2. collapse --------------------------------------------------------------
must collapse python3 "$DEF" collapse --samplesheet "$SHEET" --out "$OUT" $DROP_UTR
done_after collapse && exit 0
COLL="$OUT/proteomes"

# 3. database: extract the ~45 defensome HMMs once ------------------------
# Scanning against the subset gives identical calls (E-values depend on the
# number of target sequences, not on how many HMMs are searched, and --cut_ga
# thresholds are per family) in a small fraction of the time.
DB="$OUT/db"; SUB="$DB/defensome.hmm"
if [[ -z "$PFAM" && "$DOWNLOAD" -eq 0 && ! -s "$SUB" ]]; then
  must setup python3 "$DEF" setup --fetch --db-dir "$DB"
elif [[ "$DOWNLOAD" -eq 1 && ! -s "$SUB" ]]; then
  must setup python3 "$DEF" setup --download --db-dir "$DB"
elif [[ -n "$PFAM" ]] && { [[ ! -s "$SUB" ]] || [[ "$PFAM" -nt "$SUB" ]]; }; then
  must setup python3 "$DEF" setup --pfam "$PFAM" --db-dir "$DB"
fi

# 4. scan ------------------------------------------------------------------
if [[ -s "$SUB" ]]; then
  ARCH=()
  if [[ "$FULL_ARCH" -eq 1 ]]; then
    if [[ -n "$PFAM" ]]; then ARCH=(--full-pfam "$PFAM")
    else echo "--full-architecture needs --pfam pointing at the full Pfam-A; skipping it"; fi
  fi
  must scan python3 "$DEF" scan --proteomes "$COLL" --pfam "$SUB" --out "$OUT" \
       --threads "$THREADS" "${ARCH[@]}"
else
  echo "reusing the existing scan output in $OUT/hmmsearch"
fi
done_after scan && exit 0

# 5. core ------------------------------------------------------------------
must annotate python3 "$DEF" annotate --proteomes "$COLL" --out "$OUT"
done_after annotate && exit 0
# proteins the ORF caller discarded, recovered from transcripts when the sheet
# names any; it precedes qc so the metallothionein verdict can credit it
if python3 -c 'import csv,sys; r=list(csv.DictReader(open(sys.argv[1]),delimiter="\t")); sys.exit(0 if any((x.get("nucleotide") or "").strip() for x in r) else 1)' "$SHEET"; then
  step short-orfs python3 "$DEF" short-orfs --samplesheet "$SHEET" --out "$OUT"
fi
step qc       python3 "$DEF" qc       --out "$OUT"
step extract  python3 "$DEF" extract  --out "$OUT" --proteomes "$COLL"
step domains  python3 "$DEF" domains  --out "$OUT"

# 5. CYP clans and trees ---------------------------------------------------
if [[ -n "$CYPREFS" && -f "$CYPREFS" ]]; then
  step cyp python3 "$DEF" cyp --out "$OUT" --refs "$CYPREFS" --family CYP --threads "$THREADS"
fi
if python3 "$DEF" doctor --require trees >/dev/null 2>&1; then
step trees python3 "$DEF" trees --out "$OUT" --threads "$THREADS" \
      --families "${DEF_TREE_FAMILIES:-CYP,CCE,GST,UGT,ABC_full}"
else
  echo; echo "skipping trees: MAFFT/FastTree not available (python3 $DEF doctor)"
fi

# 6. the comparison --------------------------------------------------------
# compare only when the sheet really describes several methods per species
NMETH=$(python3 -c 'import csv,sys; r=list(csv.DictReader(open(sys.argv[1]),delimiter="\t")); print(len({(x.get("method") or "").strip() for x in r} - {""}))' "$SHEET")
if [[ "${NMETH:-0}" -ge 2 ]]; then
  step compare python3 "$DEF" compare --samplesheet "$SHEET" --out "$OUT" \
        --reference "$REF" --complete-only ${EXCLUDE:+--exclude "$EXCLUDE"}
else
  echo; echo "skipping compare: the sample sheet has $NMETH annotation method(s)."
  echo "  It pairs several annotations of the SAME species; with one proteome per"
  echo "  species there is nothing to pair. Use --metadata/--group-by instead."
fi
done_after compare && exit 0

# 7. report and dashboard --------------------------------------------------
if [[ -n "$META" && -n "$GROUP" ]]; then
  step report python3 "$DEF" report --out "$OUT" --metadata "$META" --group-by "$GROUP"
else
  step report python3 "$DEF" report --out "$OUT"
fi
step dashboard python3 "$DEF" dashboard --out "$OUT" --proteomes "$COLL" \
      ${META:+--metadata "$META"} ${GROUP:+--group-by "$GROUP"} \
      ${TREE:+--species-tree "$TREE"}

echo; echo "########## SUMMARY $(date) ##########"
if [[ ${#FAILED[@]} -eq 0 ]]; then
  echo "all steps OK"
else
  echo "FAILED: ${FAILED[*]}"
fi
echo "results:    $OUT"
echo "comparison: $OUT/comparison/"
echo "dashboard:  $OUT/dashboard/index.html"
[[ ${#FAILED[@]} -eq 0 ]] || exit 1
