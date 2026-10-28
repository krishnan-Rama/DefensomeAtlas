# External datasets

Everything below is now built into `defensome.py`. There are no helper scripts
to lose, and no relative paths to get wrong:

```bash
python3 defensome.py dataset-help          # how to obtain and use each one
python3 defensome.py cyp-refs --help       # build a clan reference
python3 defensome.py chem-triage --help    # triage P450RDB
```

## You have not downloaded ICPD yet

`db/` contains abc_membrane, abc_tran, coe, cyp, flybase, gst_*, interactions,
p450rdb and pfam. There is no `icpd/`. Get it from
http://www.insectp450.net/ui/#/page/download (Protein Sequence, Sequence
information table, and References2), put it in `db/icpd/`, then:

```bash
python3 defensome.py cyp-refs --source icpd \
    --fasta db/icpd/ICPD_protein.fasta \
    --table db/icpd/ICPD_sequence_info.xlsx \
    --orders Lepidoptera --out db/icpd/ref.faa
cd-hit -i db/icpd/ref.faa -o db/icpd/ref.c80.faa -c 0.8 -n 5
python3 defensome.py cyp --out results/ --refs db/icpd/ref.c80.faa
```

## Until then, use what you already have

```bash
python3 defensome.py cyp-refs --source flybase \
    --fasta db/flybase/flybase_dmel_clan_ref.canonical.faa \
    --out db/flybase/ref.faa
python3 defensome.py cyp --out results/ --refs db/flybase/ref.faa
```

This is a Drosophila-only reference, roughly 350 Myr from Lepidoptera, so
expect LOW-confidence clan calls. It keeps you running; it is not the answer.

## P450RDB

```bash
python3 defensome.py chem-triage --dir db/p450rdb --out db/p450rdb/triage
```

Reports kingdom and species composition, counts insect rows, melts the parallel
substrate columns into enzyme-substrate pairs while dropping O2, H2O and the
other cofactors, and writes `plant_metabolite_inventory.tsv`.

Two things to expect. The visible records are plant and microbial, and the
reactions are BIOSYNTHETIC: casbene to ketocasbene, tryptophan to an oxime in
cyanogenic glycoside biosynthesis, oleanolic acid to hederagenin. A plant CYP79
building a cyanogenic glycoside and a lepidopteran CYP6 degrading one are
opposite directions of chemistry. And a reaction database has NO NEGATIVES:
absence of a record is not evidence of a non-substrate, so no binary classifier
can be trained on it honestly.

Its real value is the plant metabolite inventory, joinable to your host-plant
table. That is the chemical challenge each moth species actually faces.

## The human dataset

Ni et al. (2025) Sci Data 12:1427 is the only source with curated
NON-substrates, which gives it two jobs, neither of them predicting insect
biology: reproduce their GCN as a methods benchmark (published MCC 0.51 to
0.72), then measure how far a human-trained model transfers to functionally
validated insect pairs. Only CYP3A4 shares a clade with insect CYP6/CYP9, at
roughly 20-25% identity, so expect a large fall. That number is the result.

Licence is CC BY-NC-ND. NoDerivatives constrains merged datasets; check the
Figshare record separately and get it in writing.

## The curation decision that matters most

ICPD's References2 links genes to phenotypes. Separating
*overexpressed in a resistant strain* (correlative, and not evidence the enzyme
metabolises anything) from *heterologously expressed and shown to metabolise
compound X* (functional) is the single most important step. Only the second is
a substrate label, and reporting how many records survive the split is a
contribution on its own.
