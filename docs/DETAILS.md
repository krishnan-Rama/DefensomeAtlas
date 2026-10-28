# defensome in detail

The [README](../README.md) gets you from proteomes to a dashboard. This page is
the reasoning behind it: what the comparisons can and cannot support, how the
calls are made, and how the dashboard behaves at catalogue scale.

## The comparison this is built for

If you have the same species annotated more than one way, the sample sheet
turns that into a measured variable:

| sample_id | species | method | fasta | gene_regex |
|---|---|---|---|---|
| A_aura_GG_040725_okay | A_aura | genome_guided | … | `_i\d+\w*$` |
| A_aura_300725_okay | A_aura | denovo | … | `_i\d+\w*$` |
| A_aura_pep_DToL | A_aura | DToL | … | *(gene: tag)* |

`compare` then pairs every contrast **within species**, so unlike a
cross-species comparison it is not confounded by phylogeny: the same genome is
being measured twice. It reports per-family recovery against a reference
method, Wilcoxon signed-rank tests with BH correction, and paired figures.

Isoform collapsing is not optional, and it is not uniform. On a real
ten-species set (EvidentialGene `okay` transcriptome sets against DToL genome
annotations), collapsing every protein to its gene gave a median of 7.6% for
de novo, 6.6% for genome-guided, and anywhere from 0.8% to 54.1% for DToL
depending on which pipeline produced the annotation. EvidentialGene's `okay`
sets have already removed most redundant isoforms, which is why the
transcriptome figures are low; an unfiltered Trinity assembly collapses far
more.

**Three further things that data showed, and that the pipeline now handles:**

- *"DToL" is not one method.* In that set it was BRAKER for four species
  (median about 18,900 genes after collapsing), Ensembl for five (about
  14,000) and raw AUGUSTUS for one. `inspect` records the `pipeline` and
  `compare` reports recovery split by it.
- *Trinity "genes" are not biological genes.* A fragmented assembly splits one
  real gene across several `c_g` components, so after collapsing isoforms the
  transcriptome arms still ran a median 1.97x above the genome annotation of the
  same species (range 1.13-2.42x). `compare --complete-only` counts only intact domain
  architectures, which strips fragments.
- *Per-10k normalisation is wrong for a paired comparison.* Within one species
  the true defensome is fixed. Dividing by an annotation-dependent gene total
  penalises exactly the arms with the noisiest totals. `compare` therefore
  defaults to raw counts; per-10k is opt-in and meant for comparing different
  species.

---


## Small proteins, and the families that came back empty

Three families returned zero in early runs. They had three different causes,
and a count table cannot tell them apart, so every zero now gets a labelled
verdict in `qc_zero_families.tsv` rather than a shrug.

**Metallothioneins are not a Pfam problem, they are a missing-protein problem.**
Insect MTs are 40–64 residues with 10–13 cysteines (Drosophila MtnA–F), and
InterPro places them in family 5, Diptera (IPR000966), which lists **no Pfam
member at all**. On top of that, TransDecoder keeps ORFs of at least 100
residues by default, so a transcriptome-derived proteome **cannot contain
one**. No domain threshold can find a protein that is not in the file.
`inspect` now reports the shortest protein in each proteome and flags
`NO_SHORT_PROTEINS` when an ORF caller has cut them off; `short-orfs` goes
back to the transcripts, takes ATG-to-stop ORFs in the 25–99 residue band the
caller discarded, and applies a composition rule.

**MATE was a threshold problem.** Its hits sit just below Pfam's gathering
threshold in divergent taxa. Families marked `rescue=yes` in the map get a
second search with an E-value cutoff instead of `--cut_ga`, and each hit must
then pass the map's `validate` rules (`copies>=2` for MATE, since a real
transporter has two MatE repeats; `len<=120` and `cys>=0.15` for MT). Rescued
calls are written to `counts_rescued.tsv` and `rescue_calls.tsv` with the
reason for every pass and fail, and are **never merged into `counts.tsv`**.

**Some families were never searched.** Their accession is absent from the Pfam
release being used. `setup` writes a `MANIFEST.tsv` recording exactly which of
the map's accessions your release contains, so this is ruled out first instead
of being mistaken for a strict map rule.

The composition screen for metallothioneins rejects on **aromatic residues**,
not cysteine alone. Defensins and knottins are also short and cysteine-rich but
carry F, W or Y; a cysteine-only rule accepted a defensin decoy while rejecting
a genuine MT at 23% cysteine. On 331,342 short ORFs from random sequence the
current rule returns zero candidates, and it recovers Drosophila MtnB exactly
from a transcript in a 50 Mb assembly in about ten seconds.

None of this is merged into the counts. It is evidence, presented beside the
counts in the dashboard's **Annotation methods** tab, and it needs confirming
by alignment before you cite a number.

## Running on SLURM

Three jobs chained by dependency: a collapse job, a scan array with one task
per proteome, and one analysis job for everything downstream.

```bash
# once, on the login node: generate the sample sheet and CHECK it
python3 defensome.py inspect --proteomes /path/to/peps --out results

# edit slurm/config.sh for your partitions, account and modules, then
bash slurm/submit.sh --samplesheet results/samplesheet.tsv \
     --pfam /db/Pfam-A.hmm --out results --drop-utrorf --dry-run
bash slurm/submit.sh --samplesheet results/samplesheet.tsv \
     --pfam /db/Pfam-A.hmm --out results --drop-utrorf
```

`slurm/config.sh` is the only file you should need to edit. It holds the
partitions, resources, and two module-loading functions. There are two because
HMMER and a modern Python stack usually need incompatible compiler toolchains
on EasyBuild/Lmod clusters, and loading both makes Lmod swap GCCcore under one
of them. `load_scan` pairs HMMER with a Python on the *same* GCCcore (the scan
needs no third-party Python packages); `load_analysis` records absolute paths
to HMMER, MAFFT and FastTree first, then loads the Python stack last.

Before submitting anything, `submit.sh` runs a preflight in *both*
environments using those same functions, and installs matplotlib to `--user`
if the analysis Python lacks it. It fails on the login node rather than
thirty jobs later.

**Resuming.** Rerun the identical command. Samples with a finished scan are
skipped, and if every collapsed proteome exists the collapse job is skipped
too. The array reads `results/scan_list.txt`, written at submission, so the
task-to-sample mapping is fixed and auditable.

Inside the analysis job, `annotate` is fail-fast: if it breaks, the job stops
with one message instead of a cascade. Independent steps such as `trees` and
`cyp` record a failure and let the rest finish.

### CYP clans are not optional

`submit.sh` refuses to proceed without a clan-labelled CYP reference, because a
missing one used to mean clans were skipped with a single log line. Pass
`--cyp-refs`, set `CYP_REFS` in `slurm/config.sh`, or pass `--no-cyp` to waive
it deliberately. To add clans to a finished run without rescanning:

```bash
bash slurm/submit.sh --samplesheet results/samplesheet.tsv --pfam /db/Pfam-A.hmm \
     --out results --cyp-refs ref.c80.faa --steps cyp,genetree,dashboard
```

`--steps` submits only the analysis job and runs only the named steps.

## The dashboard is a site

`dashboard` writes a folder, `results/dashboard/`, not a single file:

```
dashboard/
  index.html            the comparative views; opens on its own
  data/species/*.js     one per species: genes, domain coordinates, sequences
  data/families/*.js    one per family: every member gene
  data/trees/*.js       trees too large to embed, loaded when opened
  data/fasta/*.faa.gz   every called protein, one file per family
```

Open `index.html` from disk; no server is needed, because the data files are
loaded with script tags, which browsers allow from `file://` where they block
`fetch()`. Copy or host the folder as it is. `index.html` alone still opens and
every comparative view works; only species, gene and family-member sections
need `data/`, and they say so when it is missing.

It is organised like a database front end, not a report:

* **Home**: the species and families at a glance, cards for each group, and
  what each group's defensome is made of.
* **Species**: a sortable list; each species has a page with its taxonomy as a
  breadcrumb, its traits, links out (Wikipedia, GBIF, NCBI Taxonomy, GoaT, ENA),
  its defensome against the typical species, and every gene with its domains
  drawn to scale. FASTA for the species, or one family in it, is a button.
* **Families**: every family with a plain name and its InterPro entries. Each
  family has a page: distribution by group, ranked species, every member gene,
  its trees, and downloads.
* **Gene pages**: the sequence, domain architecture, the closest characterised
  proteins (BLAST against Swiss-Prot), the other copies in that genome, the
  trees that hold it, and **structure**, below.
* **Compare**: pick any trait (habitat, diet, class...) and see every family's
  median in every group as one grid; pick a group to see what differs from the
  rest. Two species side by side is one click further down.
* **Explore**: heatmap, ordination, correlations, CYP clans.
* **About**: everything about how the numbers were made, folded away:
  downloads, methods, gene-calling rules, quality control, domain completeness,
  isoform collapse, group tests.

The pages carry no explanatory text; the **?** button in the header turns it
on. **Search** (press `/`) finds species by binomial or common name, families
by name, plain name or Pfam accession, taxa at any rank, and protein IDs.

Traits come from the metadata file: any column whose values repeat (habitat,
diet, feeding guild, life stage) becomes a filter, a grouping, and a Compare
option, with no configuration. Add a `diet` column and it is there.

Every view has its own address (`#/species/<id>`, `#/family/CYP`,
`#/gene/<id>/<protein>`), so the back button works and a link can be sent.

Sequences come from `--proteomes` (read once per species) or, without it, from
the family FASTA that `extract` wrote. `--no-sequences` leaves them out. The
species data is rewritten only when an input changes; `--force` rewrites it.
Only the defensome Pfam models are searched, so a gene page shows those
domains and no others; the page says so.

## Structure and binding for a gene you choose

Every gene page has a **Structure** panel. The family's cofactor is preset
(haem for P450s, FAD and NADP+ for FMOs, glutathione for GSTs, ATP for ABC
transporters...), because a substrate pose in an empty P450 pocket means
nothing. Add a ligand as SMILES (the PubChem button looks it up by name) and
download a [Boltz-2](https://github.com/jwohlwend/boltz) input or an AlphaFold
Server input, or copy the command to run it on the cluster:

```bash
python3 defensome.py structure --out results --gene 'SAMPLE|PROTEIN' \
    --ligand 'CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1' --ligand-name permethrin --run
# without --run it writes results/structures/<gene>__<ligand>/run_boltz.sh to submit on a GPU node
```

Rerun `dashboard` and the prediction appears on the gene page: the model in 3D
coloured by confidence (pLDDT), pTM, and for a ligand Boltz-2's binder
probability and predicted IC50. The 3D viewer is fetched when opened, so that
one panel needs a connection; everything else stays offline.

Read the binding numbers for what they are: a predicted affinity from a
predicted complex. They rank ligands for one enzyme reasonably; they are not
evidence that the enzyme metabolises the compound, which needs turnover, not
binding.

## Families with domains Pfam's thresholds miss

For families the map marks `rescue` (MATE, metallothionein), copies found
below Pfam's gathering threshold that pass the family's own checks (two MatE
repeats; a short cysteine-rich protein) are real genes the threshold missed.
The dashboard counts them and marks them "below threshold" on every gene table
and gene page. `counts.tsv`, and the statistics built on it, stay strictly at
the gathering threshold; `dashboard --ga-only` makes the dashboard match it.

A zero is not always a miss. Lončar et al. (2016) found no MATE in
*Drosophila*, nematode or sea urchin genomes, so zero MATE in insects can be
real absence; fish, where MATE is present, are the positive control. The
About page's rules section gives every empty family a verdict.

## The tree view

The dashboard's Trees tab is built for reading gene trees as evidence about
annotation method. Copies of one gene recovered by DToL, genome-guided and de
novo proteomes should sit together; a clade built from one method alone is a
method-specific gene model, a fragment family, or a contaminant. Click a tip for
its record and a link to its gene page; a gene page links back to every tree
that holds it, opening the tree on that gene.

- **Rings** for method, species, CYP clan, domain completeness, and any trait,
  each in its own palette so no colour means two things. Legends carry counts.
- **Single-source clades**: every maximal clade of three or more tips sharing
  one method (or one species, or one clan), with a button that zooms to it.
- **Click a branch point** for that clade's composition by method, species,
  clan and completeness, and a flag when it is single-method. Download its tips.
- **Scroll to zoom, drag to pan**, radial or rectangular, cladogram or
  phylogram with a scale bar, FastTree support shown on demand, search that
  highlights matching tips. Labels appear once they would be legible.
- PNG and SVG export of the current view.

A 3,000-tip tree renders in about 130 ms: all branches are one path, each ring
colour is one path, and hit testing is geometric rather than per-element.

`test/test_trees.sh` checks this against a fixture with a known answer (twelve
orthologous clades plus three injected de novo-only clades that must be found
exactly), then opens the dashboard in a real headless browser, moves a real
mouse to computed tip positions before and after wheel zoom and drag pan in
both layouts, and asserts the tooltip names the right tip.



## Running across distant lineages

The map's defaults were tuned on insects. If your proteomes span phyla (fish,
amphibians, molluscs and insects in one run), three things need attention
before the numbers mean anything.

**Metallothionein needs two accessions.** `PF00131` (Metallothio) covers
vertebrate class-I MTs and should be found at the gathering threshold in fish
and amphibians. No Pfam model covers insect MTs, so the rescue pass, the
composition screen and `short-orfs` are the detectors there. Both accessions
are now in the map; confirm their names in `MANIFEST.tsv` for your release.

**MATE is a positive control in a mixed run.** MATE/SLC47 is present in fish
and orthologous to human ([zebrafish](https://doi.org/10.1038/srep28937)). If
MATE is called at the gathering threshold in your fish but only ever *rescued*
below it in your insects, that contrast is informative: the rescued insect
calls need scrutiny, and the fish calls validate the family rule. A mixed run
gives you this for free; an insect-only run does not.

**CYP clan assignment does not transfer across lineages.** Best-hit assignment
has no "unassigned" outcome, so an insect CYP2/CYP3/CYP4/MITO reference will
place a fish or mollusc P450 into an insect clan with a large and meaningless
margin. `cyp` now prints this warning and takes `--samples` to restrict the run
to the lineage the references came from.

Also remember that per-10k normalisation compares gene budgets, and a pufferfish
and a mayfly last shared an ancestor around 600 Ma. Cross-phylum defensome
differences are confounded by genome size, annotation pipeline, gene-family
turnover and ecology at once. The paired within-species design in `compare` is
the only comparison in this pipeline that is not confounded by phylogeny.

## Running at catalogue scale

Measured on 1,575 proteomes (a full Darwin Tree of Life release), 600 defensome
proteins each, single core:

| step | time | peak RSS |
|---|---|---|
| inspect | 6 s | 85 MB |
| annotate | ~190 s | ~1 GB |
| qc | 1 s | 97 MB |
| extract | 4 s | 433 MB |
| domains | 101 s | 599 MB |
| report | 191 s | 2.4 GB |
| dashboard | minutes (reads each proteome once) | ~2 GB |

`index.html` stays small (the comparative matrices, about 10 MB at 1,575
species) because nothing per gene is embedded. The per-gene data costs roughly
0.4 kB per gene with sequences, so about half a gigabyte for 1.4M gene calls,
plus the gzipped family FASTA. A rerun with unchanged inputs skips that pass.
The slowest view in the browser is the PCA, at about 4 s. `trees` is the exception: it is
not run at this scale, because a family FASTA then holds tens of thousands of
sequences and MAFFT is the limit, not this pipeline.

Give `report` at least 4 GB. Figures that size themselves per species are
clamped, so a 1,575-species heatmap is one tall block rather than a request for
a 630-inch canvas, which is what used to get the process killed by the OOM
killer.

## Filtering and grouping by taxonomy

With a metadata file the dashboard carries **every** categorical column, not
just `--group-by`. A filter bar above the tabs gives one dropdown per
taxonomic rank, and they cascade: pick `class = Insecta` and the order list
narrows to insect orders. Filtering rebuilds a view of the payload, so every
tab, figure and export shows the subset, and the header states what is
selected. A second selector regroups the whole dashboard by any rank without
regenerating the file.

Any column qualifies, not only the recognised ranks: a column is treated as a
grouping variable when its values repeat, which is what distinguishes `habitat`
from a measurement. So `habitat` and `life_stage_habitat` in a metadata file
become filters and grouping options beside `phylum` and `class`, and the
heatmap can be sorted on any of them independently of the grouping, letting you
group by class and still read the rows in habitat order.

## The Overview

The landing view is a figure, not a table: mean copies per species in each
group, split by the part of the defence system they belong to, as a share or in
absolute copies, and split by category or by family. It follows the filter bar
and the "Group by" selector, so the same panel answers "what does the Insecta
defensome look like" and "how does marine differ from freshwater" without
regenerating anything. The per-species list lives under Species, with a page
behind every name.

## Comparing two taxa

The Compare tab asks the same question at two levels. Two genomes is the
original: one species against another, family by family. Two taxa treats every
species in each group as a replicate, so each family gets a median per taxon, a
Mann-Whitney test and a Benjamini-Hochberg adjusted q. Any metadata column
works, so `habitat = marine` against `habitat = freshwater` and
`class = Insecta` against `class = Actinopteri` are the same control.

Species are not independent. A family that survives the correction is a
candidate for a phylogenetic model, not a result; the test counts 35 marine
species as 35 replicates, which they are not.

## Reading a dashboard of 1,575 species

Three things change once the species axis is longer than the page:

* **Names.** Sample ids carry the suffix of the file they came from
  (`Takifugu_rubripes_pep`). That is provenance, not a name, and it is stripped
  from every axis, legend and table. A `common_name` column, if the metadata has
  one, is shown after the binomial as *Takifugu rubripes* (torafugu), with the
  binomial italic and the vernacular upright. The tree search matches the
  vernacular too, so "zebrafish" finds *Danio rerio*.
* **Labels.** Above 60 rows the heatmap stops naming every row. If the rows are
  grouped it labels each contiguous band with the group and its size instead;
  otherwise it thins the names to a readable subset. The ordination labels the
  30 points furthest from the centroid, skipping any that would land on a label
  already placed, and says how many of the total it drew. Everything unlabelled
  is still on hover.
* **Bar charts.** The per-species bars in the family tab show the top 40 with a
  button to open all of them, rather than a 27,000-pixel scroll by default.

`common_name` is filled by `make_taxonomy_metadata.py` from the taxdump's
"genbank common name", and only for species matched exactly: a genus-level
match would put a genus-level vernacular beside a binomial, which is wrong.

## Trees at catalogue scale

A family FASTA from 1,575 proteomes holds tens of thousands of sequences and
MAFFT is the limit. `trees --per-species N` keeps the longest N sequences per
species when a family exceeds `--max-seqs`, giving a tree with the same
taxonomic span that actually builds (1,575 species x 2 = 3,150 tips rather than
~80,000). Without it such families are skipped, and the log says so.

## Trees for a question, not only for a family

A whole-family tree across 1,575 species answers few questions. `trees` takes
restrictions, names the output after them so a subset never overwrites the
whole-family tree, and records what went in beside it:

```bash
# CYPs of the Lepidoptera only
python3 defensome.py trees --out results --families CYP \
    --metadata results/taxonomy.tsv --where order=Lepidoptera
#   -> trees/CYP__order-lepidoptera.tre  (+ .json sidecar)

# marine and freshwater fish, intact gene models only
python3 defensome.py trees --out results --families GST --metadata results/taxonomy.tsv \
    --where habitat=marine,freshwater --where class=Actinopteri --complete-only

# the CYP3 clan, where insect detoxification P450s expand
python3 defensome.py trees --out results --families CYP --clan CYP3

# Lepidoptera CYPs, one representative species per family, every copy kept
python3 defensome.py trees --out results --families CYP --metadata results/taxonomy.tsv \
    --where order=Lepidoptera --thin family:1 --complete-only
```

**Thin species, not genes.** 442 Lepidoptera hold about 40,000 CYPs, far more
than MAFFT and FastTree will take. `--thin family:1` keeps one species per
lepidopteran family (the one whose copy number is closest to its family's
median, so neither an inflated nor a sparse assembly) and every copy in it: the
tree still shows the whole CYP repertoire across the order's span.
`--per-species N` does the opposite, keeping N copies per species, and a tree of
three CYPs out of ~90 per genome says little about the family; it remains as a
last resort and now keeps the copies closest to the median length.

**Over-long models are dropped** before alignment: anything longer than
`--max-len-ratio` (default 2) times the median length, which in a P450 set means
read-through or fused gene models. They open gaps across every other sequence,
and MAFFT's memory grows with the square of the longest sequence. The log and
the sidecar record how many went. If MAFFT is still killed, the run says it ran
out of memory and moves on to the next family instead of stopping.

`--where` matches any metadata column, case-insensitively; commas inside one
`--where` are OR, repeated `--where` are AND. `--complete-only` drops PARTIAL
and FRAGMENT models (from `domains`), which otherwise pull alignments apart and
sit on long branches. `--clan` needs `cyp`. `--name` overrides the generated
label. Rerun `dashboard` and every tree appears on its family's page and in the
Trees view, labelled with its restrictions.

## Grouping a catalogue-scale run by taxonomy

```bash
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir taxdump && tar xzf taxdump.tar.gz -C taxdump
# straight from the proteome directory, before anything else has run:
python3 scripts/make_taxonomy_metadata.py /path/to/proteomes taxdump results/taxonomy.tsv
# or from a sample sheet once `inspect` has written one:
python3 scripts/make_taxonomy_metadata.py results/samplesheet.tsv taxdump results/taxonomy.tsv
python3 defensome.py dashboard --out results --metadata results/taxonomy.tsv --group-by class
```

Offline, no API key, no rate limit. It writes `sample_id`, `species`,
`matched`, `common_name` and `superkingdom`..`family`. Unmatched species are
kept with blank ranks so you can see them rather than lose them. Without a grouping the species axis
of a 1,575-row heatmap is in arbitrary order and shows nothing.

## Metadata: the grouping is yours to supply

`compare` pairs annotations of one species. When you have one proteome per
species the contrast is between species, and that needs a grouping variable
nothing can infer. Build a skeleton and fill it in:

```bash
python3 scripts/make_metadata_template.py results/samplesheet.tsv results/metadata.tsv habitat class
```

`examples/metadata_35species.tsv` is a worked one for a mixed aquatic set.
Grouping is optional: without `--metadata` the pipeline runs end to end and
simply omits the group means, group tests and trait colouring.

`qc`'s control expectations are priors for insect genomes. On other taxa supply
your own with `qc --expect family<TAB>count`, or read the verdicts as advisory;
the per-species outlier test beside them is self-calibrating.

## Prerequisites: what to upload, what the pipeline fetches

Nothing needs to come off your HPC. The repository is self-contained except
for one database, and the pipeline can fetch that itself.

Per-family HMM files you may already have (one directory per family, each
`hmmpress`ed) are **not** used by this pipeline and do not need to be: `setup`
builds a single multi-model file covering all ~45 map accessions, which
`hmmsearch` searches in one pass. A directory holding six of them cannot
replace it.

UniProt proteomes work too, and are often the only source for laboratory
model species: `inspect` reads their `GN=` gene tags, detects the `sp|`/`tr|`
accession format, and reports how many entries are unreviewed. Collapsing
matters more there than anywhere else, since a UniProt proteome bundles
isoforms (zebrafish UP000000437 is 68,855 proteins for roughly 26,000 genes).
Mixing UniProt with Ensembl and BRAKER in one run adds a third annotation
source, so record it as a covariate. UniProt's proteome table also ships a
BUSCO completeness figure per proteome, which is a better annotation-quality
covariate than anything this pipeline can compute, and worth copying into your
metadata.

| what | size | how |
|---|---|---|
| Pfam-A.hmm | ~300 MB gz | only to rebuild the subset or for `--full-architecture`: `--pfam /path/to/Pfam-A.hmm[.gz]` or `--download-pfam` |
| the defensome HMM subset | ~1 MB | downloaded prebuilt from the GitHub release by `run`/`setup`, or built by `setup --pfam`; **this** is what gets searched |
| CYP clan reference | ~1 MB | build from ICPD with `cyp-refs`; see [DATASETS.md](DATASETS.md) |
| HMMER, MAFFT, FastTree | – | `envs/environment.yml`, or cluster modules |

**Do not commit Pfam-A to the repository.** `setup` extracts only the ~45
accessions the map uses into `<out>/db/defensome.hmm`, about 1 MB, and writes
`MANIFEST.tsv` recording the release. Searching the subset gives identical
calls, because `hmmsearch` scores each HMM independently and `--cut_ga` uses
per-family thresholds, while E-values depend on the number of target
sequences rather than on how many models are searched. It is also dramatically
faster: the full-Pfam architecture pass is now optional (`--full-architecture`)
and runs only on proteins that already hit a defensome family.

The one file worth committing is your clan reference if you want other people
to reproduce your CYP clans exactly. It is about 1 MB and gzips well. Commit it
under `db/` with the ICPD citation (Wu et al. 2025, *Mol Ecol Resour* 25:e14070)
in the README, and note ICPD's terms.

## Subcommands

| | |
|---|---|
| `doctor` | report which tools and libraries are present and which steps can run |
| `setup` | extract the defensome HMMs from Pfam-A and record the release |
| `short-orfs` | recover proteins the ORF caller discarded, from the nucleotide FASTA in the sample sheet |
| `inspect` | examine headers, propose a gene-ID rule, write a draft sample sheet |
| `collapse` | longest protein per gene, per the sample sheet |
| `scan` | `hmmsearch` against Pfam-A with `--cut_ga` |
| `annotate` | apply the family architecture rules, build count matrices |
| `qc` | annotation quality from low-copy control families |
| `compare` | paired comparison across annotation methods (raw by default; `--complete-only`, `--exclude`, `--per10k`) |
| `domains` | per-protein domain architecture and completeness |
| `extract` / `trees` | per-family FASTA; MAFFT + FastTree per family or on a subset (`--where`, `--complete-only`, `--clan`) |
| `cyp` / `cyp-refs` / `cyp-benchmark` | CYP clan assignment, reference building, hold-out accuracy |
| `icpd-evidence` / `chem-triage` | triage external functional datasets |
| `report` / `dashboard` | figures, and the dashboard site with species, family and gene pages |
| `assets` / `dataset-help` | write embedded files out, dataset guidance |

---

## The family map

`defensome_map.tsv`, one row per family, editable:

| column | meaning |
|---|---|
| `pfam_ids` | comma-separated Pfam accessions |
| `rule` | `ANY` = at least one; `ALL` = every one, in the same protein |
| `min_cov` | minimum fraction of the HMM aligned, over merged non-overlapping segments |
| `min_len` | minimum protein length |
| `tier` | `CORE` = safe to compare; `BROAD` = huge superfamily, counts inflated |

`ALL` is what stops a lone GST_N domain being called a GST, or a bare P-loop
NTPase being called an ABC transporter. `python3 defensome.py assets --write .`
writes the map out; a file next to the script overrides the embedded copy.

---

## Design decisions worth knowing

**Raw counts versus per-10k.** `counts.tsv` holds the number of genes matching
each family. `counts_per10k.tsv` holds that number divided by the total genes
in the same proteome, times 10,000: a share, not a count. It removes exactly
one nuisance variable, proteome size, and nothing else. It is not z-scoring
(the heatmap and PCA do that separately, purely for colour) and it is not a
statistical normalisation in the RNA-seq sense.

Which one is correct depends on the comparison:

- **Different species** → per-10k. Raw counts partly track genome size: across
  a 43-species set the Spearman correlation between raw CORE defensome and
  proteome size was 0.51, and −0.37 after dividing.
- **Same species, different annotation methods** → raw. The true defensome is
  fixed, so dividing by an annotation-dependent denominator injects the
  artefact you are measuring. `compare` therefore defaults to raw.

And per-10k cannot repair a bad denominator. A genome annotation left at 70,619
genes makes its own defensome look four times smaller; normalisation propagates
that rather than fixing it. Check `qc` and the `inspect` gene counts first.

The denominator is genes once `collapse` has run and proteins otherwise, and
every label now says which, because those are not the same number when a
proteome carries several isoforms per gene.

**Low-copy families are a free annotation-quality score.** `qc` compares each
family's median against a literature expectation (a map problem) and each
species against the observed median (a genome problem). A genome with seven
times the median catalase count has duplicated gene models. This replaces a
BUSCO duplication run for most purposes.

**Families with no variance are reported, not hidden.** An all-zero family and
a family constant at exactly 1 are both findings. They get their own panel.

**Sequences over 100,000 aa are excluded before searching.** HMMER's pipeline
aborts with SIGABRT on any target that long rather than skipping it, so one
merged-scaffold gene model in one proteome kills the whole scan. Genome
annotations do emit them. The excluded proteins are listed in
`oversized_proteins.tsv`, never silently dropped.

**HMM coverage is summed over merged non-overlapping segments.** Models that
match in two pieces (FMO's paired Rossmann folds, ABC_membrane) would otherwise
be scored on their best segment and mislabelled truncated.

---

## Testing

```bash
bash test/test_pipeline.sh               # needs no external tools
bash test/test_site.sh                   # node + playwright for the browser part
bash test/test_dashboard.sh results/dashboard/index.html   # needs node
```

`test_pipeline.sh` builds a proteome set with a **known answer**: exactly 120
genes per transcriptome sample over a variable number of isoforms, 150 per
genome sample, CYPs injected at 20/16/12 so recovery is exactly 0.80 and 0.60,
and 8 real GST pairs plus 6 unpaired GST_N decoys the `ALL` rule must reject.
It then asserts every one of those numbers.

`test_dashboard.sh` extracts the payload and app code back out of the generated
HTML, runs them under a minimal DOM shim in node, executes every view, and
checks the maths: z-scores have mean 0 and sd 1, the PCA eigensolver is
verified against a matrix with a known answer, and the Newick parser is checked
for branch lengths, pruning, and not mistaking a support value for a tip name.

`test_site.sh` builds 140 species with every gene called, decoys in each
proteome, and two-domain families. It checks that a `--where` tree holds exactly
the genes of that taxon, that AND, OR and `--complete-only` combine correctly,
and that an unknown column stops the run. It then opens the site from `file://`
in a real browser: search, taxon and breadcrumb navigation, the back button,
lazily loaded genes, members and trees, and that every sequence a gene page or
species FASTA shows is byte-identical to the proteome and excludes the decoys.

---

## Citation

See `CITATION.cff`. Please also cite HMMER, Pfam, MAFFT, FastTree and
Goldstone et al. (2006) for the defensome concept.

## Licence

MIT. See `LICENSE`.
