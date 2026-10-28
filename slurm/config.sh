#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# Cluster configuration. This is the ONLY file you should need to edit.
# Everything else in slurm/ sources it.
# ---------------------------------------------------------------------------

# Partitions. The scan array is embarrassingly parallel and short per task;
# the analysis job is one long serial job.
PART_SCAN="${PART_SCAN:-epyc}"
PART_ANALYSIS="${PART_ANALYSIS:-epyc_long}"

# Account. Leave empty if your cluster does not use one: an empty
# --account= is rejected the same way a wrong one is, so it is omitted.
# Find yours with: sacctmgr show assoc user=$USER format=Account -P
ACCOUNT="${ACCOUNT:-}"

# Clan-labelled CYP reference for `cyp`. Without one the pipeline cannot
# assign CYP clans, and submit.sh will refuse to proceed silently. Build one
# from ICPD with `defensome.py cyp-refs` (see docs/DATASETS.md), or point at
# any FASTA whose headers carry clan=CYP2|CYP3|CYP4|MITO.
CYP_REFS="${CYP_REFS:-}"

# Resources
SCAN_CPUS="${SCAN_CPUS:-16}"
SCAN_MEM="${SCAN_MEM:-16G}"
SCAN_TIME="${SCAN_TIME:-12:00:00}"

# Largest array a single submission may have. Clusters cap queued jobs per
# partition (Skomer: epyc 500, epyc_long 20), and an array above the cap is
# rejected outright. submit.sh picks a chunk size to stay under this.
MAX_ARRAY="${MAX_ARRAY:-400}"
SCAN_CONCURRENT="${SCAN_CONCURRENT:-12}"      # max array tasks at once
ANALYSIS_CPUS="${ANALYSIS_CPUS:-16}"
ANALYSIS_MEM="${ANALYSIS_MEM:-64G}"
ANALYSIS_TIME="${ANALYSIS_TIME:-48:00:00}"

# ---------------------------------------------------------------------------
# Module loading, one function per job type.
#
# The two functions exist because HMMER and a modern Python stack usually
# need incompatible compiler toolchains on EasyBuild/Lmod clusters. Loading
# them together makes Lmod silently swap GCCcore underneath one of them.
#
#   load_scan      HMMER plus a Python built on the SAME GCCcore. `scan`
#                  needs only the Python standard library, no pandas.
#   load_analysis  records absolute paths to HMMER/MAFFT/FastTree first,
#                  then loads the Python stack last. defensome.py calls the
#                  recorded binaries directly (DEFENSOME_* variables), so the
#                  toolchain swap no longer matters.
#
# The defaults below are for an EasyBuild cluster with these modules:
#   HMMER/3.4-gompi-2023b           (GCCcore 13.2.0)
#   Python/3.11.5-GCCcore-13.2.0    (GCCcore 13.2.0, matches HMMER)
#   MAFFT/7.526-GCC-14.2.0-with-extensions, FastTree/2.2-GCCcore-14.2.0
#   CD-HIT/4.8.1-GCC-14.3.0, SciPy-bundle/2025.07-gfbf-2025b
# Run `module avail` and change the version strings for your cluster.
# ---------------------------------------------------------------------------
_quiet_load() { module load "$@" 2>/dev/null || module load "$@"; }

load_scan() {
  module purge
  _quiet_load HMMER/3.4-gompi-2023b
  _quiet_load Python/3.11.5-GCCcore-13.2.0
}

load_analysis() {
  module purge
  _quiet_load HMMER/3.4-gompi-2023b
  export DEFENSOME_HMMSEARCH="$(command -v hmmsearch)"
  export DEFENSOME_HMMBUILD="$(command -v hmmbuild)"
  _quiet_load MAFFT/7.526-GCC-14.2.0-with-extensions FastTree/2.2-GCCcore-14.2.0
  export DEFENSOME_MAFFT="$(command -v mafft)"
  export DEFENSOME_FASTTREE="$(command -v FastTree)"
  _quiet_load CD-HIT/4.8.1-GCC-14.3.0
  export DEFENSOME_CD_HIT="$(command -v cd-hit)"
  # Python stack LAST, so whatever it swaps cannot affect the paths above
  _quiet_load SciPy-bundle/2025.07-gfbf-2025b
}

# Clusters without Lmod/EasyBuild: replace both functions with, for example,
#   load_scan()     { source ~/miniforge3/bin/activate defensome; }
#   load_analysis() { source ~/miniforge3/bin/activate defensome; }
# using envs/environment.yml, which has every tool in one environment.
