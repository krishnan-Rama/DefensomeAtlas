# DEFENSOME

**Find and compare the chemical defensome (the genes animals and plants use to
detoxify, transport and survive chemicals) across any set of proteomes, and
browse the result as an interactive website.**

Give it a folder of protein FASTA files, one per species. It finds every
cytochrome P450, GST, UGT, ABC transporter, carboxylesterase, metallothionein
and ~30 other defence gene families, counts them per species, checks each gene
model, and builds a dashboard where you can search species, compare groups
(habitat, diet, taxonomy), open any gene to see its sequence and domains, and
set up a structure prediction with a ligand.

[![CI](https://github.com/krishnan-Rama/DEFENSOME/actions/workflows/ci.yml/badge.svg)](https://github.com/krishnan-Rama/DEFENSOME/actions)

---

## Contents

1. [Install](#1-install)
2. [Run it](#2-run-it)
3. [Choose your species](#3-choose-your-species)
4. [Add taxonomy and traits](#4-add-taxonomy-and-traits)
5. [Open the dashboard](#5-open-the-dashboard)
6. [Optional extras: trees, CYP clans, structures](#6-optional-extras)
7. [Running on a cluster](#7-running-on-a-cluster-slurm)
8. [What comes out](#8-what-comes-out)
9. [The HMM database](#9-the-hmm-database)
10. [Troubleshooting](#10-troubleshooting)

---

## 1. Install

### Get the code

```bash
git clone https://github.com/krishnan-Rama/DEFENSOME.git
cd DEFENSOME
python3 defensome.py --version
```

There is nothing to compile and nothing to `pip install` for the tool itself:
`defensome.py` is one file, with the family map and the dashboard built in.

### Get the dependencies (pick one route)

| | needed for | |
|---|---|---|
| Python ≥ 3.9 with **pandas, numpy** | everything | required |
| **HMMER 3** (`hmmsearch`) | finding the genes | required |
| matplotlib, scipy | figures, group tests | recommended |
| MAFFT, FastTree | gene trees, CYP clans | optional |
| Boltz-2 (GPU) | structure / ligand prediction | optional, separate environment |

**A. conda or mamba** (laptop, workstation or cluster; no admin rights needed)

```bash
conda env create -f envs/environment.yml
conda activate defensome
```

**B. Cluster modules** (module names differ between clusters; `module avail` lists yours)

```bash
module load HMMER MAFFT FastTree Python
pip install --user -r requirements.txt
```

**C. macOS with Homebrew**

```bash
brew install hmmer mafft fasttree
pip3 install -r requirements.txt
```

Then check what can run:

```bash
python3 defensome.py doctor
```

### Optional: Boltz-2 for structures

Boltz needs its own environment: it pins older numpy and scipy than the main
tool, so use Python 3.10–3.12 (on newer Pythons pip tries to compile scipy and
fails). It also needs a GPU to be practical.

```bash
conda env create -f envs/boltz.yml            # or:
python3.11 -m venv ~/envs/boltz && ~/envs/boltz/bin/pip install boltz
```

You never activate it for `defensome.py`; you only put it on the `PATH` of the
GPU job that runs a prediction (see [structures](#structures-and-ligand-binding)).

---

## 2. Run it

```bash
python3 defensome.py run --proteomes my_proteomes/ --out results
```

That is the whole pipeline: it collapses isoforms to one protein per gene,
downloads the defensome HMMs (a few MB, once), searches every proteome, calls
and checks the genes, and writes the dashboard to
`results/dashboard/index.html`.

**Your proteomes:** one protein FASTA per species (`.faa`, `.fa` or `.fasta`,
gzipped or not). The file name becomes the species ID, so name files by
species: `Danio_rerio.faa`, `Bombyx_mori_pep.fa.gz`.

**Useful options**

| option | what it does |
|---|---|
| `--species LIST` | run only these species ([below](#3-choose-your-species)) |
| `--threads 8` | CPUs for the search (default 4) |
| `--metadata FILE --group-by COL` | your taxonomy or traits ([below](#4-add-taxonomy-and-traits)) |
| `--taxdump DIR` | build the taxonomy automatically from NCBI |
| `--trees CYP,GST` | also build gene trees for these families |
| `--cyp-refs FILE` | assign P450s to clans (CYP2, CYP3, CYP4, mito) |
| `--no-collapse` | your files already have one protein per gene |
| `--review` | stop after the sample sheet, to check how genes were grouped |
| `--hmm FILE` / `--pfam FILE` | use your own HMMs instead of downloading them |
| `--no-figures` | skip the static PNG/PDF figures |

Run `python3 defensome.py run --help` for the full list. A run can be repeated:
finished searches are skipped, so adding species or changing options only
redoes what changed.

---

## 3. Choose your species

Point `--proteomes` at a folder holding many species and pick some with
`--species`, as a comma list or a file with one name per line:

```bash
python3 defensome.py run --proteomes all_proteomes/ --out moths \
    --species "Bombyx mori,Manduca sexta,Danaus plexippus"

python3 defensome.py run --proteomes all_proteomes/ --out fish --species fish.txt
```

```
# fish.txt: one per line; lines starting with # are ignored
Danio rerio
Takifugu_rubripes_pep
Salmo                    # a genus name selects every Salmo species
```

A name matches a file by its full name, or as a binomial or genus followed by
the rest of the file name, ignoring case, spaces and underscores. A name that
matches nothing stops the run before any work, and says which one.

---

## 4. Add taxonomy and traits

Groups are what make the comparison: *which families are bigger in marine
species, in herbivores, in Lepidoptera?*

**Automatic taxonomy.** Download the NCBI taxonomy once (about 60 MB; no
account needed) and pass it in:

```bash
wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mkdir taxdump && tar xzf taxdump.tar.gz -C taxdump
python3 defensome.py run --proteomes my_proteomes/ --out results --taxdump taxdump/
```

This adds phylum, class, order and family, plus common names
(*Danio rerio* (zebrafish)), and groups by class.

**Your own traits.** Write a tab-separated file with one row per proteome and
a column holding the file names (without extension). Every other column whose
values repeat becomes a filter, a grouping and a comparison in the dashboard:

```
sample_id            class        order         habitat      diet
Danio_rerio          Actinopteri  Cypriniformes freshwater   omnivore
Bombyx_mori_pep      Insecta      Lepidoptera   terrestrial  herbivore
Crassostrea_gigas    Bivalvia     Ostreida      marine       filter
```

```bash
python3 defensome.py run --proteomes my_proteomes/ --out results \
    --metadata traits.tsv --group-by habitat
```

`scripts/make_metadata_template.py` writes a blank template from a sample
sheet; `scripts/make_taxonomy_metadata.py` is the taxdump step on its own.

---

## 5. Open the dashboard

Open `results/dashboard/index.html` in a browser. No server, no internet
(except for the 3D structure viewer).

- **Home**: your species and families at a glance, and what each group's defensome is made of.
- **Species**: every species; each has a page with its taxonomy, links out (Wikipedia, GBIF, NCBI, GoaT, ENA), its defensome against the typical species, and every gene.
- **Families**: every gene family; each has a page with its spread across groups, members, trees and downloads.
- **Gene pages**: sequence (copy or download FASTA), domains drawn to scale, BLAST links, and structure prediction.
- **Compare**: pick a trait, see every family's median in every group, then one group against the rest.
- **Explore**: heatmap, ordination, correlations, CYP clans.
- **About**: downloads, methods text, gene-calling rules and quality control.

Press `/` to search; the **?** button shows explanations.

**Keep the `data/` folder beside `index.html`**: species, gene and family pages
load from it. To share, zip the whole `dashboard/` folder. If you work on a
cluster, copy it to your computer first:

```bash
rsync -av user@cluster:/path/to/results/dashboard/ ./dashboard/
```

Keep it out of OneDrive or Dropbox while browsing: it is thousands of small
files and sync clients lag behind.

---

## 6. Optional extras

### Gene trees

```bash
python3 defensome.py run ... --trees CYP,GST            # during the run, or afterwards:
python3 defensome.py trees --out results --families CYP --complete-only
python3 defensome.py trees --out results --families CYP --metadata traits.tsv \
    --where order=Lepidoptera --thin family:1           # one species per moth family
python3 defensome.py trees --out results --families CYP --clan CYP3
```

`--where` restricts by any metadata column, `--complete-only` drops fragments,
`--thin family:1` keeps one representative species per family with all its
copies. Large families need memory: run trees as a job with `--mem=32G`.
Rerun `dashboard` afterwards and trees appear on the family page.

### CYP clans

```bash
python3 defensome.py dataset-help                       # where to get reference sequences
python3 defensome.py run ... --cyp-refs cyp_clan_refs.faa
```

### Structures and ligand binding

Every gene page has a **Structure** panel: the family's cofactor is preset
(haem for P450s, FAD for FMOs, glutathione for GSTs, ATP for ABC transporters),
you add a ligand as SMILES, and download an input for
[Boltz-2](https://github.com/jwohlwend/boltz) or
[AlphaFold Server](https://alphafoldserver.com), or copy the command. From the
command line:

```bash
# write the input (on any machine)
python3 defensome.py structure --out results --gene 'Bombyx_mori_pep|BMSK0001234' \
    --ligand 'CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1' --ligand-name permethrin

# predict on a GPU node, with the Boltz environment on PATH
sbatch -p gpu --gres=gpu:1 --mem=32G --time=2:00:00 --wrap \
  "export PATH=$HOME/envs/boltz/bin:\$PATH; bash results/structures/<folder>/run_boltz.sh"

# show it on the gene page
python3 defensome.py dashboard --out results --proteomes results/proteomes
```

`--gene` takes `SAMPLE|PROTEIN`; the easiest way to get it exactly right is
the **command** button on that gene's page. The first run downloads the Boltz
weights (several GB, to `~/.boltz`) and every run fetches an alignment from
the ColabFold server, so the GPU node needs internet.

Predicted binding is a ranking aid, not evidence of metabolism. Check that the
cofactor and ligand sit in the active site before trusting the numbers.

---

## 7. Running on a cluster (SLURM)

For dozens to thousands of proteomes, `slurm/submit.sh` runs the search as a
job array and everything else as one job:

```bash
python3 defensome.py inspect --proteomes /path/to/proteomes --out results   # once; check results/samplesheet.tsv
# edit slurm/config.sh: partitions, account, module names
bash slurm/submit.sh --samplesheet results/samplesheet.tsv --out results \
     --metadata results/taxonomy.tsv --group-by class --no-cyp --dry-run
bash slurm/submit.sh --samplesheet results/samplesheet.tsv --out results \
     --metadata results/taxonomy.tsv --group-by class --no-cyp
```

The HMMs are fetched on the login node before submission. Rerun the same
command to resume; `--steps dashboard` reruns one step. See
[`slurm/README.md`](slurm/README.md).

---

## 8. What comes out

```
results/
  dashboard/index.html      the website (keep data/ beside it)
  counts.tsv                copies per family per species
  counts_per10k.tsv         the same per 10,000 genes
  gene_calls_all.tsv        every gene: species, family, protein, length, coverage
  domains/                  domain coordinates and completeness of every gene
  fasta/<family>.faa        protein sequences per family
  qc_annotation.tsv         annotation quality per species
  figures/                  PNG/PDF figures and report_all_figures.pdf
  trees/                    gene trees and what went into each
  structures/               Boltz inputs and predictions
  samplesheet.tsv           how each proteome was read
  db/                       the HMMs used, and the Pfam release they came from
```

`counts.tsv` uses only genes at Pfam's curated gathering thresholds. The
dashboard also counts, and marks, validated genes just below threshold for
short or divergent families (MATE, metallothionein); `dashboard --ga-only`
turns that off.

---

## 9. The HMM database

DEFENSOME searches ~45 Pfam models, not all 20,000. They are downloaded
automatically from this repository's
[`pfam-db` release](https://github.com/krishnan-Rama/DEFENSOME/releases/tag/pfam-db)
as `defensome_db.tar.gz`, checked against its SHA-256, and saved with the Pfam
release number in `results/db/PFAM_VERSION`.

| situation | do this |
|---|---|
| normal use | nothing; `run` fetches it |
| no internet on the machine | download `defensome_db.tar.gz` elsewhere, then `--db-url /path/to/defensome_db.tar.gz`, or put it in `DEFENSOME/db/` |
| a newer Pfam, or your own copy | `--pfam /path/to/Pfam-A.hmm.gz` (builds the subset in a few minutes) |
| you edited the family map | rebuild with `--pfam`; `run` warns if the map uses a model the database lacks |

The family map is `defensome_map.tsv`: one row per family with its Pfam
accessions, the domain rule, coverage and length thresholds. Edit it to add a
family; `python3 defensome.py assets` writes a fresh copy out. Maintainers
publishing a new database: [`docs/RELEASING.md`](docs/RELEASING.md).

---

## 10. Troubleshooting

| message | cause and fix |
|---|---|
| `hmmsearch not found` | load or install HMMER; `python3 defensome.py doctor` lists what is missing |
| `no proteome matches X` | check the spelling against the file names it lists |
| `could not download the defensome HMMs` | offline; see [section 9](#9-the-hmm-database) |
| MAFFT `Killed` in trees | out of memory: run as a job with more `--mem`, or use `--thin`, `--clan`, `--where` |
| `rsync: connection unexpectedly closed` with `stty` warnings | your cluster `~/.bashrc` runs `stty`; change it to `[[ $- == *i* ]] && stty ...` |
| Boltz install builds scipy and fails | Python is too new; make the Boltz environment with Python 3.11 |
| Boltz: `No module named 'cuequivariance_torch'` | already handled (`--no_kernels`); for an old `run_boltz.sh`, add `--no_kernels` |
| a family is zero everywhere | About → Gene-calling rules gives each empty family a verdict; a zero can be real absence |

**More:** the reasoning, scale behaviour and design decisions are in
[`docs/DETAILS.md`](docs/DETAILS.md); external CYP datasets in
[`docs/DATASETS.md`](docs/DATASETS.md); version changes in
[`UPGRADING.md`](UPGRADING.md). Every step is also its own subcommand
(`python3 defensome.py --help`).

## Citation

See [`CITATION.cff`](CITATION.cff). Please also cite HMMER, Pfam/InterPro,
MAFFT, FastTree, and Boltz-2 if you use structures.

## Licence

MIT. The Pfam models in the HMM database are CC0 (public domain).
