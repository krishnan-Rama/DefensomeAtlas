#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Submit the whole pipeline to SLURM. Three jobs, chained by dependency:
#
#   00_collapse   one job      isoform collapse, per the sample sheet
#   01_scan       array job    one hmmsearch per collapsed proteome
#   02_analysis   one job      annotate -> qc -> domains -> cyp -> trees ->
#                              compare -> report -> dashboard
#
# Resubmitting resumes: samples already scanned are skipped, and if every
# collapsed proteome already exists the collapse job is skipped too.
#
# Usage:
#   bash slurm/submit.sh --samplesheet results/samplesheet.tsv \
#        --pfam /db/Pfam-A.hmm --out results [options]
#
# First time? Generate the sample sheet on the login node and CHECK it:
#   python3 defensome.py inspect --proteomes /path/to/peps --out results
# ---------------------------------------------------------------------------
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
source "$HERE/slurm/config.sh"

SHEET=""; PFAM=""; OUT=""; EXCLUDE=""; CYPREFS=""; TREE=""; META=""; GROUP=""
REF="DToL"; DROP=""; DRY=0; FORCE_COLLAPSE=0; STEPS=""; NO_CYP=0; DOWNLOAD=0; FULL_ARCH=0
CHUNK=0; MAX_ARRAY="${MAX_ARRAY:-400}"

usage() { sed -n '2,20p' "$0"; cat <<'U'
Options:
  --samplesheet FILE   required
  --pfam FILE          Pfam-A.hmm or .hmm.gz (or --download-pfam)
  --download-pfam      fetch Pfam-A from EBI into <out>/db instead
  --full-architecture  also record non-defensome domains (needs the full --pfam)
  --out DIR            required
  --exclude LIST       species or sample_ids to leave out of compare (e.g. A_lych)
  --reference METHOD   recovery reference method (default DToL)
  --drop-utrorf        discard EvidentialGene utrorf entries at collapse
  --cyp-refs FILE      clan-labelled CYP reference FASTA
  --species-tree FILE  Newick species tree
  --metadata FILE --group-by COL
  --concurrent N       max simultaneous scan tasks (default from config.sh)
  --chunk N            proteomes per array task (default: chosen so the array
                       stays under MAX_ARRAY, currently 400)
  --force-collapse     rerun collapse even if proteomes exist
  --steps LIST         rerun only these analysis steps, skipping collapse and
                       scan: annotate,qc,extract,domains,
                       cyp,trees,genetree,compare,report,tree,dashboard
                       (e.g. --steps trees,dashboard)
  --no-cyp             proceed without CYP clan assignment (not recommended)
  --dry-run            print what would be submitted, submit nothing
Edit slurm/config.sh for partitions, account, resources and modules.
U
exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SHEET=$2; shift 2;;
    --pfam) PFAM=$2; shift 2;;
    --out) OUT=$2; shift 2;;
    --exclude) EXCLUDE=$2; shift 2;;
    --reference) REF=$2; shift 2;;
    --drop-utrorf) DROP="--drop-utrorf"; shift;;
    --cyp-refs) CYPREFS=$2; shift 2;;
    --species-tree) TREE=$2; shift 2;;
    --metadata) META=$2; shift 2;;
    --group-by) GROUP=$2; shift 2;;
    --concurrent) SCAN_CONCURRENT=$2; shift 2;;
    --chunk) CHUNK=$2; shift 2;;
    --force-collapse) FORCE_COLLAPSE=1; shift;;
    --download-pfam) DOWNLOAD=1; shift;;
    --full-architecture) FULL_ARCH=1; shift;;
    --steps) STEPS=$2; shift 2;;
    --no-cyp) NO_CYP=1; shift;;
    --dry-run) DRY=1; shift;;
    -h|--help) usage;;
    *) echo "unknown option: $1" >&2; usage;;
  esac
done

fail() { echo "ERROR: $*" >&2; exit 1; }
[[ -n "$SHEET" && -n "$OUT" ]] || usage
[[ -f "$SHEET" ]] || fail "sample sheet not found: $SHEET"
if [[ -n "$PFAM" ]]; then [[ -f "$PFAM" ]] || fail "Pfam file not found: $PFAM"; fi
if [[ -n "$META" ]]; then [[ -f "$META" ]] || fail "metadata not found: $META"; fi
# A missing reference used to mean clans were skipped with one log line that
# nobody reads. Now it is a hard stop unless explicitly waived.
: "${CYPREFS:=$CYP_REFS}"
if [[ -n "$CYPREFS" ]]; then
  [[ -f "$CYPREFS" ]] || fail "CYP refs not found: $CYPREFS"
elif [[ "$NO_CYP" -eq 0 ]]; then
  cat >&2 <<'NOCYP'
ERROR: no CYP clan reference given, so CYP clans would not be assigned.

  Pass --cyp-refs /path/to/ref.c80.faa, or set CYP_REFS in slurm/config.sh.
  Build one from ICPD with:
    python3 defensome.py cyp-refs --source icpd --fasta all_prot.fa \
            --table information.xlsx --orders Lepidoptera --out ref.faa
    cd-hit -i ref.faa -o ref.c80.faa -c 0.8 -n 5
  Or rerun with --no-cyp to proceed without clans deliberately.
NOCYP
  exit 1
fi
if [[ -n "$TREE" ]]; then [[ -f "$TREE" ]] || fail "species tree not found: $TREE"; fi

# absolute paths: jobs do not reliably inherit the working directory
mkdir -p "$OUT"; OUT=$(cd "$OUT" && pwd); mkdir -p "$OUT/logs"
# `readlink -f` yields an EMPTY STRING when a parent directory is missing, so
# a mistyped or placeholder path would silently become "argument not given".
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
SHEET=$(abs_file "--samplesheet" "$SHEET") || exit 1
PFAM=$(abs_file "--pfam" "$PFAM") || exit 1
META=$(abs_file "--metadata" "$META") || exit 1
CYPREFS=$(abs_file "--cyp-refs" "$CYPREFS") || exit 1
TREE=$(abs_file "--species-tree" "$TREE") || exit 1
DB="$OUT/db"
# Neither --pfam nor --download-pfam: use the prebuilt defensome HMMs. Fetch
# them here, on the login node, because compute nodes are often offline.
if [[ -z "$PFAM" && "$DOWNLOAD" -eq 0 && ! -s "$DB/defensome.hmm" && "$DRY" -eq 0 ]]; then
  python3 "$HERE/defensome.py" setup --fetch --db-dir "$DB" \
    || fail "could not fetch the prebuilt HMMs; pass --pfam /path/to/Pfam-A.hmm.gz instead"
fi

# --- preflight, on THIS node, with the SAME module functions the jobs use ---
echo "########## preflight ##########"
( load_scan >/dev/null 2>&1
  command -v hmmsearch >/dev/null || { echo "scan env: hmmsearch NOT found"; exit 1; }
  command -v python3 >/dev/null   || { echo "scan env: python3 NOT found"; exit 1; }
  python3 "$HERE/defensome.py" --version >/dev/null || { echo "scan env: defensome.py will not run"; exit 1; }
  echo "  scan env      hmmsearch $(command -v hmmsearch)"
  echo "                python3   $(python3 -c 'import sys;print(sys.version.split()[0])')"
) || fail "load_scan() in slurm/config.sh does not give a working environment"
( load_analysis >/dev/null 2>&1
  python3 -c 'import pandas, numpy' 2>/dev/null || { echo "analysis env: pandas/numpy missing"; exit 1; }
  echo "  analysis env  python3   $(python3 -c 'import sys;print(sys.version.split()[0])') with pandas"
  python3 -c 'import matplotlib' 2>/dev/null || {
    echo "  analysis env: matplotlib missing, installing to --user (figures need it)"
    python3 -m pip install --user --quiet matplotlib openpyxl scipy || \
      echo "  WARNING: install failed; figures will be skipped"; }
) || fail "load_analysis() in slurm/config.sh does not give pandas and numpy"

# --- the sample list the array will index into ---
mapfile -t SAMPLES < <(python3 - "$SHEET" <<'PY'
import csv, sys
for r in csv.DictReader(open(sys.argv[1]), delimiter="\t"):
    print(r["sample_id"])
PY
)
[[ ${#SAMPLES[@]} -gt 0 ]] || fail "no samples in $SHEET"

LIST="$OUT/scan_list.txt"; : > "$LIST"
DONE=0
for s in "${SAMPLES[@]}"; do
  if [[ -s "$OUT/hmmsearch/$s.domtblout.gz" ]]; then DONE=$((DONE+1)); else echo "$s" >> "$LIST"; fi
done
N=$(wc -l < "$LIST")

# A cluster limits how many jobs may be queued at once (Skomer: epyc 500,
# epyc_long 20). A one-task-per-proteome array of 1,575 is simply rejected, so
# each task takes a chunk of the list instead.
if [[ "$CHUNK" -le 0 ]]; then
  CHUNK=$(( (N + MAX_ARRAY - 1) / MAX_ARRAY ))
  [[ "$CHUNK" -lt 1 ]] && CHUNK=1
fi
NTASK=$(( (N + CHUNK - 1) / CHUNK ))

NEED_COLLAPSE=$FORCE_COLLAPSE
for s in "${SAMPLES[@]}"; do
  [[ -s "$OUT/proteomes/$s.faa" ]] || NEED_COLLAPSE=1
done

ACC=(); [[ -n "$ACCOUNT" ]] && ACC=(--account="$ACCOUNT")
EXPORTS="ALL,DEF_HOME=$HERE,DEF_OUT=$OUT,DEF_SHEET=$SHEET,DEF_PFAM=$PFAM"
EXPORTS+=",DEF_DROP_UTR=$DROP,DEF_EXCLUDE=$EXCLUDE,DEF_REFERENCE=$REF"
EXPORTS+=",DEF_CYP_REFS=$CYPREFS,DEF_SPECIES_TREE=$TREE,DEF_META=$META,DEF_GROUP=$GROUP"
EXPORTS+=",DEF_STEPS=${STEPS//,/:},DEF_DB=$DB,DEF_DOWNLOAD=$DOWNLOAD,DEF_FULL_ARCH=$FULL_ARCH"
EXPORTS+=",DEF_CHUNK=$CHUNK"

cat <<SUMMARY

sample sheet  $SHEET  (${#SAMPLES[@]} samples)
pfam          $PFAM
out           $OUT
hmm database  $DB/defensome.hmm $([[ -s "$DB/defensome.hmm" ]] && echo "(exists, setup skipped)" || echo "(will be built from ${PFAM:-$([[ $DOWNLOAD -eq 1 ]] && echo "an EBI download" || echo "the prebuilt release")})")
cyp refs      ${CYPREFS:-<none: clans will NOT be assigned>}
steps         ${STEPS:-all}
collapse      $([[ -n "$STEPS" ]] && echo "skipped: --steps given" || { [[ $NEED_COLLAPSE -eq 1 ]] && echo "will run" || echo "skipped: every proteome already collapsed"; })
scan          $N to scan, $DONE already done
              array of $NTASK task(s) x $CHUNK proteome(s), %$SCAN_CONCURRENT at a time, on $PART_SCAN
analysis      on $PART_ANALYSIS, ${ANALYSIS_CPUS} cpus, ${ANALYSIS_TIME}
exclude       ${EXCLUDE:-<none>}
account       ${ACCOUNT:-<none, flag omitted>}
SUMMARY

if [[ "$DRY" -eq 1 ]]; then
  echo; echo "--dry-run: nothing submitted. Scan list:"; cat "$LIST"; exit 0
fi

# --steps: rerun part of the analysis only; nothing upstream is resubmitted
if [[ -n "$STEPS" ]]; then
  [[ -s "$OUT/counts.tsv" || ",$STEPS," == *",annotate,"* ]] || \
    fail "--steps without annotate needs an existing $OUT/counts.tsv"
  J2=$(sbatch --parsable "${ACC[@]}" --partition="$PART_ANALYSIS" \
       --cpus-per-task="$ANALYSIS_CPUS" --mem="$ANALYSIS_MEM" --time="$ANALYSIS_TIME" \
       --output="$OUT/logs/analysis_%j.out" --export="$EXPORTS" \
       "$HERE/slurm/02_analysis.sh")
  echo; echo "submitted analysis  $J2  (steps: $STEPS)"
  echo "logs:    $OUT/logs/analysis_$J2.out"
  exit 0
fi

# collapse must finish before scan starts, and scan must know the array size
# now; the list comes from the sheet, so it does not depend on collapse output.
DEP=""
if [[ ! -s "$DB/defensome.hmm" ]]; then
  JS=$(sbatch --parsable "${ACC[@]}" --partition="$PART_ANALYSIS" \
       --output="$OUT/logs/setup_%j.out" --export="$EXPORTS" \
       "$HERE/slurm/00_setup.sh")
  DEP="--dependency=afterok:$JS"
  echo; echo "submitted setup     $JS"
fi

if [[ $NEED_COLLAPSE -eq 1 ]]; then
  J0=$(sbatch --parsable "${ACC[@]}" --partition="$PART_ANALYSIS" $DEP \
       --output="$OUT/logs/collapse_%j.out" --export="$EXPORTS" \
       "$HERE/slurm/00_collapse.sh")
  DEP="--dependency=afterok:$J0"
  echo; echo "submitted collapse  $J0"
fi

if [[ "$N" -gt 0 ]]; then
  J1=$(sbatch --parsable "${ACC[@]}" --partition="$PART_SCAN" $DEP \
       --array=1-"$NTASK"%"$SCAN_CONCURRENT" --cpus-per-task="$SCAN_CPUS" \
       --mem="$SCAN_MEM" --time="$SCAN_TIME" \
       --output="$OUT/logs/scan_%A_%a.out" --export="$EXPORTS" \
       "$HERE/slurm/01_scan.sh")
  # afterany, not afterok: in an array of hundreds of tasks a single failure
  # would otherwise leave the analysis job stuck on DependencyNeverSatisfied
  # for ever. The analysis job checks completeness itself and says what is
  # missing, which is far more useful than a job that silently never runs.
  DEP="--dependency=afterany:$J1"
  echo "submitted scan      $J1  ($NTASK tasks, $CHUNK proteome(s) each)"
fi

J2=$(sbatch --parsable "${ACC[@]}" --partition="$PART_ANALYSIS" $DEP \
     --cpus-per-task="$ANALYSIS_CPUS" --mem="$ANALYSIS_MEM" --time="$ANALYSIS_TIME" \
     --output="$OUT/logs/analysis_%j.out" --export="$EXPORTS" \
     "$HERE/slurm/02_analysis.sh")
echo "submitted analysis  $J2"
echo
echo "watch:   squeue -u \$USER"
echo "logs:    $OUT/logs/"
echo "resume:  rerun this exact command; finished samples are skipped"
