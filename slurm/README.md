# SLURM

| file | what it is |
|---|---|
| `config.sh` | **edit this**: partitions, account, resources, module loading |
| `submit.sh` | the one command; preflight, resume logic, dependency chain |
| `00_collapse.sh` | isoform collapse, one job |
| `01_scan.sh` | hmmsearch, array job, one task per proteome |
| `02_analysis.sh` | annotate through dashboard, one job |

```bash
bash slurm/submit.sh --samplesheet results/samplesheet.tsv --out results --no-cyp
```

Without `--pfam`, the prebuilt defensome HMMs are downloaded on the login node
before anything is submitted (compute nodes are often offline). Give
`--pfam /path/to/Pfam-A.hmm.gz` to build them from your own Pfam copy instead.

Each array task scans a **chunk** of proteomes, not one. Clusters cap how many
jobs may be queued (Skomer: epyc 500, epyc_long 20), so a 1,575-task array is
rejected outright. `submit.sh` picks a chunk size keeping the array under
`MAX_ARRAY` (400 by default): 1,575 proteomes become 394 tasks of 4. Override
with `--chunk N`. Resume works inside a chunk, so a task re-run skips the
proteomes it already finished.

The analysis job depends on the scan array with **`afterany`, not `afterok`**.
In an array of hundreds of tasks a single failure under `afterok` leaves the
analysis stuck on `DependencyNeverSatisfied` for ever. Instead the analysis job
checks for itself that every sample in the sheet has a non-empty scan output,
and refuses with the list of missing ones. `DEF_ALLOW_PARTIAL=1` overrides it.

Add `--dry-run` first to see the preflight and the scan list without
submitting. Rerun the same command to resume.
