#!/bin/bash
#SBATCH --job-name=def_setup
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
set -uo pipefail
source "$DEF_HOME/slurm/config.sh"
load_analysis
echo "=== setup $(date) on $(hostname) ==="
if [[ "${DEF_DOWNLOAD:-0}" == "1" ]]; then
  python3 "$DEF_HOME/defensome.py" setup --download --db-dir "$DEF_DB"
elif [[ -n "${DEF_PFAM:-}" ]]; then
  python3 "$DEF_HOME/defensome.py" setup --pfam "$DEF_PFAM" --db-dir "$DEF_DB"
else
  python3 "$DEF_HOME/defensome.py" setup --fetch --db-dir "$DEF_DB"
fi
