#!/bin/bash
#SBATCH --job-name=def_scan
set -uo pipefail
source "$DEF_HOME/slurm/config.sh"
load_scan

# submit.sh wrote scan_list.txt, one sample_id per line, before this array
# existed. Task N reads line N, so the index-to-sample mapping is fixed at
# submission and auditable afterwards.
# Each array task handles a CHUNK of the list, not one line. A cluster caps how
# many jobs may be submitted at once (epyc 500, epyc_long 20 on Skomer), so a
# 1,575-task array is simply rejected. Chunking keeps the array small.
CHUNK="${DEF_CHUNK:-1}"
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK + 1 ))
LAST=$(( FIRST + CHUNK - 1 ))
mapfile -t SIDS < <(sed -n "${FIRST},${LAST}p" "$DEF_OUT/scan_list.txt")
[[ ${#SIDS[@]} -gt 0 ]] || { echo "no samples at rows $FIRST-$LAST"; exit 0; }

echo "=== task $SLURM_ARRAY_TASK_ID  rows $FIRST-$LAST  ${#SIDS[@]} proteomes  $(hostname)  $(date) ==="
echo "python3: $(command -v python3)   hmmsearch: $(command -v hmmsearch)"
DB="${DEF_DB}/defensome.hmm"
[[ -s "$DB" ]] || { echo "ERROR: $DB missing; the setup job did not finish"; exit 1; }
ARCH=()
if [[ "${DEF_FULL_ARCH:-0}" == "1" && -n "${DEF_PFAM:-}" ]]; then ARCH=(--full-pfam "$DEF_PFAM"); fi

fail=0
for SID in "${SIDS[@]}"; do
  if [[ -s "$DEF_OUT/hmmsearch/$SID.domtblout.gz" ]]; then
    echo "-- $SID already scanned, skipping"; continue
  fi
  FAA="$DEF_OUT/proteomes/$SID.faa"
  if [[ ! -s "$FAA" ]]; then
    echo "-- ERROR $SID: collapsed proteome missing or empty: $FAA"; fail=1; continue
  fi
  t0=$SECONDS
  if python3 "$DEF_HOME/defensome.py" scan --proteomes "$DEF_OUT/proteomes" \
       --pfam "$DB" --out "$DEF_OUT" \
       --threads "${SLURM_CPUS_PER_TASK:-8}" --species "$SID" "${ARCH[@]}"; then
    echo "-- $SID ok in $((SECONDS-t0))s"
  else
    echo "-- ERROR $SID: scan failed"; fail=1
  fi
done
echo "=== task $SLURM_ARRAY_TASK_ID done, rc=$fail, $((SECONDS))s total $(date) ==="
exit $fail
