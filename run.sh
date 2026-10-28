#!/usr/bin/env bash
# Shortest possible use:  bash run.sh <proteome_dir> <outdir> [threads] [metadata.tsv group_column]
# The same as:  python3 defensome.py run --proteomes DIR --out DIR ...
set -euo pipefail
[[ $# -ge 2 ]] || { sed -n '2,3p' "$0"; exit 1; }
HERE=$(cd "$(dirname "$0")" && pwd)
ARGS=(--proteomes "$1" --out "$2" --threads "${3:-4}")
if [[ $# -ge 5 ]]; then ARGS+=(--metadata "$4" --group-by "$5"); fi
exec python3 "$HERE/defensome.py" run "${ARGS[@]}"
