#!/bin/bash
#SBATCH --job-name=def_collapse
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
set -uo pipefail
source "$DEF_HOME/slurm/config.sh"
load_analysis
echo "=== collapse $(date) on $(hostname) ==="
python3 "$DEF_HOME/defensome.py" --version
python3 "$DEF_HOME/defensome.py" collapse --samplesheet "$DEF_SHEET" --out "$DEF_OUT" $DEF_DROP_UTR
