# Upgrading from an earlier copy

Earlier versions needed `dashboard.html`, `dashboard.js` and
`defensome_map.tsv` next to `defensome.py`. If you extracted a tarball and also
have loose copies of those files in a parent directory, the loose ones shadow
the new ones and you get errors like:

```
ERROR: template not found: .../tool/dashboard.html
```

From 1.2.0 everything is embedded in `defensome.py`, so the clean fix is:

```bash
cd /path/to/tool
python3 defensome.py --version        # confirm which file you are running
```

If that prints a path you did not expect, or a version below 1.2.0, delete the
stale loose copies and keep only the extracted directory:

```bash
rm -f dashboard.html dashboard.js defensome_map.tsv defensome.py \
      01_scan.sh 02_report.sh 03_deep.sh submit.sh test_dashboard.sh
tar xzf defensome.tar.gz
cd defensome
python3 defensome.py --version
```

Nothing in `results/` is affected. Everything downstream of `hmmsearch` is
recomputed from the cached `domtblout` files in seconds.

---

# Troubleshooting an implausible gene count

`inspect` flags `IMPLAUSIBLE_GENE_COUNT` when a genome annotation still has more
than about 35,000 genes after collapsing. An insect genome should land roughly
between 10,000 and 30,000.

The usual cause is an unfiltered ab initio AUGUSTUS set. Its IDs look like
`AUGUSTUSALYP00000000001.t1`: the `.t1` suffix only groups isoforms when the
part before it is a *gene* ID (`g123.t1`, `g123.t2`). Here it is a unique
protein serial, so nothing groups and every prediction counts as a gene,
including many short spurious ORFs.

Diagnose it:

```bash
F=peps/A_lych_pep_DToL.fasta
grep -c '^>' $F                                   # total entries
head -3 $F                                        # any gene: tag at all?
awk '/^>/{if(l)print l;l=0;next}{l+=length}END{print l}' $F \
  | sort -n | awk '{a[NR]=$1}END{print "median length",a[int(NR/2)]; \
    n=0; for(i in a) if(a[i]<100) n++; print n, "proteins under 100 aa"}'
```

A large spike under 100 aa confirms spurious ORFs. Options, in order of
preference: obtain the filtered annotation from the DToL portal or Ensembl;
restrict to proteins with a BUSCO or Pfam hit; or leave the sample out of the
paired comparison with `compare --exclude A_lych`. Do not simply raise the
threshold. The family filters (`min_len`, `min_cov`) keep most spurious ORFs out
of the defensome counts, but they cannot fix a gene total inflated threefold,
and that total is what any normalisation divides by.

## 3.14.0: common names and decluttered axes

The dashboard is a presentation change only; `counts.tsv` and everything
upstream are untouched, so you do not rerun `scan` or `annotate`. To pick up
the common names you do need to rebuild the metadata file, because earlier
copies have no `common_name` column:

```bash
python3 scripts/make_taxonomy_metadata.py results/samplesheet.tsv taxdump results/taxonomy.tsv
python3 defensome.py dashboard --out results --metadata results/taxonomy.tsv --group-by class
```

One behavioural fix comes with it. The metadata reader used to admit a
non-rank column only when pandas reported its dtype as `object`. Pandas 3 reads
text columns as the `str` extension dtype, so on a recent pandas every
non-rank column was silently dropped and `habitat` never appeared as a filter
or grouping option. Selection is now on cardinality, which is the real
criterion. If your earlier dashboard showed only `phylum`, `class` and `order`
in the filter bar when your metadata had more, this is why; regenerate it with
the command above and the rest will be there.

## 3.15.0: fewer tabs, less empty ink

Presentation only again; nothing upstream of `counts.tsv` reruns. Regenerate
with the same `dashboard` command.

What changed and why:

* **Tabs with nothing in them are no longer in the nav.** A run without trees,
  domains, clans, group tests or a method comparison used to show five pages
  whose only content was a sentence explaining their own emptiness. The Methods
  tab already lists which outputs are missing and what to run for each, so that
  information is still one click away. A typical catalogue run goes from
  fifteen tabs to seven.
* **Related views were merged.** Ordination and Correlations are one tab,
  Patterns. Compare is a section of the Species tab. The defensome map is a
  section of Methods.
* **Panels no longer nest inside panels.** An inner panel drew a second border
  one inset further in, so nothing lined up down the page; inner sections are
  now separated by a rule at the same left edge.
* **Table columns that are empty in every row are dropped**, as are KPI tiles
  whose value is a dash. A column of dashes says a step was not run, which is
  what the Methods tab is for.
* **The heatmap column headers stay put while the rows scroll**, the matrix is
  capped at 620px of height rather than running to 27,000px, the label gutter
  and the header block are sized to the labels actually drawn, and the cells
  widen to use the panel. The PNG export is still the whole matrix, headers
  included.
* **Ordination labels no longer overlap or run off the panel.** Each label is
  measured, placed only where nothing is already drawn, and flipped to the left
  of its point near the right edge.
* **Where a sample sheet is present, labels use the species name from it**, with
  the method appended only when several proteomes of one species would
  otherwise read identically (`Takifugu rubripes · denovo`).

If you preferred a tab that has now disappeared, it disappeared because this
run has no data for it; give the run that output and the tab comes back.

## 3.16.0: the Overview is a figure again, and the tabs are back

Presentation only. Regenerate with the same `dashboard` command.

* **The separate tabs are restored.** Ordination, Correlations, Compare and Map
  are their own tabs again, as they were before 3.15. Tabs with no data in the
  run are still left out of the nav, and the Methods tab still says what is
  missing and what to run for it.
* **The Overview leads with a figure, not a table.** It is now the composition
  of each group's defensome: mean copies per species, split by defence category
  or by family, shown as a share or in absolute copies. It follows the filter
  bar and the "Group by" selector. The 140-row species table moved to the Data
  tab.
* **Category axes adapt to their labels.** The box-plot axis in the Families tab
  measured its longest label, lays the labels at an angle when they will not fit
  upright, widens the gutter so the leftmost one is not cut off, and thins them
  when even angled labels would collide. The correlation matrix sizes its
  gutters from the family names for the same reason.

The axis fix matters most when you group by `order` rather than `class`:
`Malacostracaformes (n=8)` next to `Mammaliaformes (n=10)` is where the old
fixed axis broke.

## 3.17.0: isoform collapse at catalogue scale, and comparing two taxa

* **"Genes after isoform collapse" was empty.** It was not empty, it was
  11,200 pixels wide at 140 species and would have been 126,000 at 1,575: a bar
  per species against a grouped-bar axis built for three species and three
  methods, with only one method in a catalogue run so nothing to group by. The
  panel now checks what the run can actually show. With one annotation method
  it charts the distribution of collapse percentage, marks the median, and says
  how many proteomes collapse far more than the rest; with several methods and
  a manageable number of species it still draws the bars. The tab is called
  "Isoform collapse" rather than "Annotation methods" when there is only one
  method, because there is no method comparison to be had.
* **Compare has a two-taxa mode**, and it comes first, because picking two
  species out of 1,575 is the specialist case. Choose any metadata column, two
  values of it, and every family gets a median per taxon, a Mann-Whitney test
  and a Benjamini-Hochberg adjusted q across families. Bars point toward
  whichever taxon has more; families surviving the correction are coloured, the
  rest are grey.

A caveat that is in the panel and belongs here too: species are not
independent, so a family that survives the correction is a candidate to test
with a phylogenetic model, not a result. The test treats 35 marine species as
35 replicates, which they are not.

## 4.0.0: the transcript layer

New, and nothing earlier changes: a run without transcripts produces exactly
what 3.17 did.

* `inspect --transcripts DIR` pairs each proteome with a transcript FASTA by
  species name, writes a `transcripts` column into the sample sheet, and
  previews linkage on three proteomes so you know what kind of transcript file
  you have before the long run.
* `transcripts` links every gene to its CDS (verified by translation), records
  model integrity, codon usage, nearest-paralog Ka/Ks and near-identical
  copies. `transcripts-summary` combines them across species.
* New outputs: `transcripts_samples.tsv`, `transcripts_family.tsv`,
  `counts_dedup.tsv`, `ks_pairs.tsv.gz`, and `transcripts/<sample>/` per
  proteome.
* The dashboard gains a Transcripts tab when those outputs exist.
* `run_pipeline.sh --transcripts DIR`; on SLURM the analysis job runs the layer
  automatically when the sample sheet names transcripts, and `--steps
  transcripts,dashboard` reruns just that part.

Requires MAFFT for Ka/Ks and near-identical duplicates (it is already loaded by
`load_analysis()` in `slurm/config.sh`). Linkage, integrity and codon usage
need only Python.

An existing results directory can gain the layer without rescanning: rebuild
the sample sheet with `inspect --transcripts`, then
`bash slurm/submit.sh ... --steps transcripts,dashboard`.

## 4.0.1: protein and transcript files from different releases

On the real DToL files, Abramis brama's proteins linked 59 of 300 to its
transcripts, all by sequence, none by ID: Ensembl cdna from a later release of
the annotation than the pep file. The genes that still match are the ones
unchanged between releases, a subset chosen by annotation history, and 4.0.0
would have analysed such a subset whenever it exceeded 20%.

A proteome is now used only when at least 80% of its genes link (a matched
pair links about 99%). Between 20% and 80%, or a sequence-only partial match,
it is reported as `release_mismatch`; below that as `not_the_source`. Neither
gets codon-level results, and the dashboard and `transcripts_samples.tsv` say
which and why. The fix is in the data: take the pep and cdna files from the
same release.

## 4.1.0: less on the screen at catalogue scale

Presentation only; regenerate the dashboard, nothing reruns.

* Box plots with more than 8 groups become one horizontal row per group,
  sorted by median, with no rotated labels. Groups with fewer than 3 species
  are left out and named in a line. Individual species are drawn only up to
  300 points, and above 7 groups every box is one colour, so unrelated groups
  no longer share a hue.
* Families tab: the group figure is full width; its numbers, the per-species
  completeness table and the domain architectures sit behind closed toggles;
  domain completeness is one bar and one line; the per-species ranking shows 20.
* Transcripts tab: one duplication measure at a time (selector), the 12 least
  complete families, the 8 highest and 8 lowest for codon bias, the rest behind
  toggles. At 1,236 species and 32 classes it went from over 6,000 pixels tall
  to about 3,100.

## 4.2.0: the count audit

The transcripts' main use here is an audit of the counts, so it is now a step
and the first thing on the Transcripts tab.

* `count-audit --out results --metadata taxonomy.tsv --group-by class` compares
  the groups three ways on the same species: raw counts, counts with
  near-identical copies collapsed, and counts of complete gene models only.
  Kruskal-Wallis per family, Benjamini-Hochberg within each version. Every
  family gets one verdict: survives; lost after deduplication; lost when only
  complete models count; lost under both; appears only after correction; or
  no difference. Writes `count_audit/audit_summary.tsv` (one row per family,
  with the highest and lowest group under each version), `audit_tests.tsv`
  and `audit_group_medians.tsv`.
* The Transcripts tab opens with the same audit, computed live, so it follows
  the filter bar and Group by. Families with no difference fold away. Ks is
  framed as a within-lineage clock; codon bias is folded away as supplementary;
  the coverage panel appears only when some species in view lack results.
* The pipeline runs it after `transcripts-summary` whenever a metadata file and
  grouping are given.

Kruskal-Wallis is implemented in the tool and matches SciPy to 1e-13; the
dashboard's version and the Python step give identical verdicts on the test
fixture.

## 5.0.0: the dashboard becomes a site; the transcript layer is gone

**Removed.** `transcripts`, `transcripts-summary`, `count-audit`, `inspect
--transcripts`, `run_pipeline.sh --transcripts`, the Transcripts tab and their
outputs. Transcript FASTA without read counts cannot give expression, and what
it could give (codon-level checks on 312 of 1,575 proteomes) did not justify a
layer of its own. Existing `transcripts_*.tsv` and `counts_dedup.tsv` files are
ignored, not deleted. `short-orfs` is unaffected: it still reads the
`nucleotide` column of the sample sheet.

**Changed.** `dashboard` writes `results/dashboard/index.html` plus
`results/dashboard/data/`, and removes the old `results/dashboard.html`.
`--dashboard-out` now names a folder. `--light` and `--full` are gone: nothing
per gene is embedded any more, so there is nothing to trim. New:
`--proteomes` (sequences for the gene pages), `--no-sequences`, `--force`.

**New.** Species, family and gene pages with sequences, domain architecture,
FASTA downloads and outbound links (NCBI Taxonomy, GBIF, InterPro, BLASTP);
global search; addressable views. `trees --where COL=VALUE --metadata FILE`,
`--complete-only`, `--clan`, `--name`, each tree with a `.json` sidecar.

**To move an existing run** (no rescan, nothing upstream changes):

```bash
python3 defensome.py dashboard --out results --metadata results/taxonomy.tsv \
    --group-by class --samplesheet results/samplesheet.tsv --proteomes results/proteomes
```

On SLURM: `bash slurm/submit.sh ... --steps dashboard`.

## 5.0.1: trees that MAFFT can actually align

* `trees --thin COLUMN[:N]` keeps N representative species per value of a
  metadata column, with every copy, instead of a few copies per species.
* Sequences longer than `--max-len-ratio` (default 2) times the median length
  are dropped before alignment; the log and the `.json` sidecar count them.
* `--per-species` keeps the copies closest to the median length, not the
  longest: the longest are the fused and read-through models.
* A MAFFT killed for memory is reported as such, and the remaining families
  still run.

## 5.1.0: a plainer dashboard, structures, traits

* The pages are browsing and pictures. QC flags, gene-call totals, isoform
  collapse, group tests, rules and downloads moved to **About**, folded away.
  Explanatory text is hidden; the **?** button shows it.
* **Compare** is top level: any trait as a grid of family medians by group,
  then one group against another or against all others.
* **Structure** on every gene page, with family cofactor presets, Boltz-2 and
  AlphaFold Server inputs, and the new `structure` command; finished Boltz runs
  appear on the gene page with confidence and predicted binding.
* Validated below-threshold genes of `rescue` families (MATE, metallothionein)
  are counted in the dashboard and marked; `dashboard --ga-only` restores
  GA-only counts. `counts.tsv` is unchanged.
* Species pages link to Wikipedia, GoaT and ENA; gene pages to a BLAST against
  Swiss-Prot for the closest characterised proteins.

## 5.1.1

`structure` now passes `--no_kernels` to Boltz, because `pip install boltz`
does not include the `cuequivariance_torch` kernels and the run otherwise stops
at the first layer. `--kernels` turns them back on once `boltz[cuda]` is
installed.

## 5.2.0: one command, any set of species, no Pfam download

* `defensome.py run --proteomes DIR --out DIR` runs everything to the
  dashboard. `--species` takes names, binomials or genera (a list or a file);
  `--taxdump` builds the taxonomy and common names; gzipped proteomes are read
  as they are. `all` is now the same command.
* The defensome HMMs are downloaded prebuilt (a few MB) from the repository's
  `pfam-db` release when no `--pfam` is given, in `run`, `setup --fetch`,
  `run_pipeline.sh` and `slurm/submit.sh`. Offline: `--db-url FILE`.
  Maintainers build the asset with `setup --pfam ... --bundle`
  (docs/RELEASING.md).
* The README is now a how-to; the reasoning moved to docs/DETAILS.md, dataset
  notes to docs/DATASETS.md, and the metadata helpers from test/ to scripts/.
* `requirements.txt`, `envs/environment.yml` (Python 3.11) and `envs/boltz.yml`.
