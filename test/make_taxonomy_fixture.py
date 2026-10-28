#!/usr/bin/env python3
"""Known-answer fixture for the taxonomic presentation of the dashboard.

140 species, every sample id carrying the "_pep" provenance suffix, across
4 phyla / 7 classes / 4 habitats, with a common name on exactly 100 of them.
That is past the 60-row limit at which the heatmap stops naming every row, so
it exercises the band labels, the label thinning and the ordination label cap,
none of which fire on the 3-species fixtures.

Counts are injected so the known answers are exact:
  * CYP is the only family that separates habitats: 60 / 40 / 20 / 10 by
    habitat in that order, with no noise, so sorting by habitat must produce
    contiguous, monotonically decreasing CYP blocks.
  * Catalase is 1 in every species, so it must be reported as invariant.
"""
import os, random, sys

PHYLA = {
    "Chordata":   {"Actinopteri": 40, "Amphibia": 14, "Mammalia": 10},
    "Arthropoda": {"Insecta": 46, "Malacostraca": 8},
    "Mollusca":   {"Bivalvia": 12},
    "Annelida":   {"Polychaeta": 10},
}
HABITATS = ["marine", "freshwater", "terrestrial", "estuarine"]
CYP_BY_HABITAT = {"marine": 60, "freshwater": 40, "terrestrial": 20, "estuarine": 10}
# The real map's full family list, so the figures are laid out against the
# number of columns a real run actually produces, not a convenient eight.
FAMS = ["CYP", "CCE", "FMO", "ALDH", "AKR", "SDR", "EH", "Cyt_b5", "GST",
        "GST_MAPEG", "UGT", "SULT", "ABC_full", "ABC_NBD", "MATE", "MFS",
        "Catalase", "SOD_CuZn", "SOD_Fe", "Peroxiredoxin", "GPx", "GCL",
        "Thioredoxin", "Glutaredoxin", "TrxR_GR", "Hsp70", "Hsp90", "Hsp20",
        "DnaJ", "Metallothionein", "Ferritin", "CDF_ZnT", "ZIP",
        "BTB_BACK_Kelch", "CncC_bZIP", "NR_xeno", "bHLH_PAS"]


# --catalogue: 32 classes of very unequal size, like the DToL set, for the
# layout tests that only bite at that scale
CATALOGUE = {
    "Arthropoda": {"Insecta": 520, "Arachnida": 30, "Malacostraca": 18, "Collembola": 4,
                   "Diplopoda": 3, "Chilopoda": 2, "Copepoda": 2},
    "Chordata": {"Actinopteri": 160, "Aves": 70, "Mammalia": 60, "Amphibia": 25,
                 "Lepidosauria": 12, "Chondrichthyes": 8, "Ascidiacea": 3, "Leptocardii": 1},
    "Mollusca": {"Gastropoda": 45, "Bivalvia": 30, "Cephalopoda": 6, "Polyplacophora": 2},
    "Annelida": {"Polychaeta": 25, "Clitellata": 14},
    "Cnidaria": {"Anthozoa": 20, "Scyphozoa": 6, "Hydrozoa": 5},
    "Echinodermata": {"Echinoidea": 8, "Asteroidea": 7, "Holothuroidea": 3, "Ophiuroidea": 3},
    "Streptophyta": {"Magnoliopsida": 120, "Liliopsida": 15, "Polypodiopsida": 4},
    "Porifera": {"Demospongiae": 5},
}


def build(root, catalogue=False, genes=False):
    global PHYLA
    if catalogue:
        PHYLA = CATALOGUE
    random.seed(11)
    os.makedirs(root, exist_ok=True)
    rows, k = [], 0
    for phylum, classes in PHYLA.items():
        for cls, n in classes.items():
            for _ in range(n):
                k += 1
                sp = f"Genus{cls[:3].lower()}{k:03d}_species{k:03d}"
                hab = HABITATS[k % 4]
                rows.append({
                    "sample_id": f"{sp}_pep",
                    "species": sp,
                    # exactly 100 carry a vernacular name; the rest are blank,
                    # which is what the taxdump path produces for most inverts
                    "common_name": f"{cls.lower()} number {k}" if k <= 100 else "",
                    "phylum": phylum, "class": cls, "order": f"{cls}formes",
                    "habitat": hab,
                })
    assert catalogue or len(rows) == 140, len(rows)

    with open(os.path.join(root, "metadata.tsv"), "w") as f:
        cols = ["sample_id", "species", "common_name", "phylum", "class", "order", "habitat"]
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r[c] for c in cols) + "\n")

    with open(os.path.join(root, "counts.tsv"), "w") as f:
        f.write("species\t" + "\t".join(FAMS) + "\n")
        for r in rows:
            v = []
            for fam in FAMS:
                if fam == "CYP":
                    v.append(CYP_BY_HABITAT[r["habitat"]])
                elif fam == "Catalase":
                    v.append(1)
                else:
                    v.append(random.randint(2, 40))
            r["counts"] = dict(zip(FAMS, v))
            f.write(r["sample_id"] + "\t" + "\t".join(map(str, v)) + "\n")


    # A catalogue run collapses isoforms but has only ONE annotation method,
    # which is the case that used to render a 126,000-pixel empty bar chart.
    with open(os.path.join(root, "collapse_stats.tsv"), "w") as f:
        f.write("sample_id\tspecies\tmethod\tn_proteins_in\tn_genes_out\t"
                "dropped_utrorf\tcollapse_pct\n")
        for r in rows:
            n_in = random.randint(14000, 62000)
            pct = round(random.uniform(0, 42), 1)
            f.write(f"{r['sample_id']}\t{r['species']}\tDToL\t{n_in}\t"
                    f"{int(n_in * (1 - pct / 100))}\t0\t{pct}\n")

    with open(os.path.join(root, "proteome_sizes.tsv"), "w") as f:
        f.write("species\tn_proteins\n")
        for r in rows:
            f.write(f"{r['sample_id']}\t{random.randint(12000, 48000)}\n")

    if genes:
        write_genes(root, rows)
    print(f"{root}: {len(rows)} species, {len(PHYLA)} phyla, {len(HABITATS)} habitats, "
          f"100 with a common name" + (", with gene-level data" if genes else ""))


AA = "ACDEFGHIKLMNPQRSTVWY"
# Two-domain families get both domains, so architectures have more than one box.
TWO_DOMAIN = {"GST": [("GST_N", "PF02798", 76), ("GST_C", "PF00043", 95)],
              "ABC_full": [("ABC_membrane", "PF00664", 274), ("ABC_tran", "PF00005", 137)]}


def seq_for(pid, n):
    """Deterministic sequence per protein, so a test can recompute it."""
    r = random.Random(pid)
    return "M" + "".join(r.choice(AA) for _ in range(n - 1))


def write_genes(root, rows):
    """Gene-level outputs consistent with counts.tsv: one gene call per copy,
    a domain table, completeness, proteomes holding those genes plus decoys,
    and two CYP trees (whole family, and one order) with their sidecars.
    Known answers:
      * every gene call's sequence is in its proteome, under its protein ID;
      * protein IDs carry "." and "-", and tree tips are sample|protein;
      * every 7th gene is PARTIAL and every 11th (not 7th) FRAGMENT."""
    r = random.Random(23)
    os.makedirs(os.path.join(root, "proteomes"), exist_ok=True)
    os.makedirs(os.path.join(root, "domains"), exist_ok=True)
    os.makedirs(os.path.join(root, "trees"), exist_ok=True)
    calls = open(os.path.join(root, "gene_calls_all.tsv"), "w")
    calls.write("species\tcategory\tfamily\tprotein\tprot_len\thmm_cov\n")
    dom = open(os.path.join(root, "domains", "domain_table.tsv"), "w")
    dom.write("species\tprotein\tprot_len\tpfam_name\tpfam_id\thmm_len\thmm_from\thmm_to\t"
              "ali_from\tali_to\tscore\tdom_cov\n")
    comp = open(os.path.join(root, "domains", "completeness.tsv"), "w")
    comp.write("species\tprotein\tfamily\tcategory\tstatus\tn_required\tn_complete\t"
               "min_domain_cov\tmean_domain_cov\tprot_len\tlength_ratio\tdomain_coverages\n")
    cyp_tips = []
    k = 0
    for row in rows:
        sid = row["sample_id"]
        tag = row["species"].split("_")[0][-6:]
        with open(os.path.join(root, "proteomes", sid + ".faa"), "w") as fa:
            for d in range(5):                       # decoys: not defensome genes
                pid = f"{tag}.decoy{d}.t1"
                fa.write(f">{pid} decoy\n{seq_for(sid + pid, 120)}\n")
            for fam in FAMS:
                for c in range(row["counts"][fam]):
                    k += 1
                    pid = f"{tag}.{fam.lower()}-g{c + 1}.t1"
                    n = r.randint(60, 260)
                    fa.write(f">{pid} gene={fam}{c}\n{seq_for(sid + pid, n)}\n")
                    calls.write(f"{sid}\tX\t{fam}\t{pid}\t{n}\t0.9\n")
                    st = "PARTIAL" if k % 7 == 0 else "FRAGMENT" if k % 11 == 0 else "COMPLETE"
                    comp.write(f"{sid}\t{pid}\t{fam}\tX\t{st}\t1\t1\t0.9\t0.9\t{n}\t1.0\t\n")
                    doms = TWO_DOMAIN.get(fam, [(fam + "_dom", "PF%05d" % (100 + FAMS.index(fam)), 120)])
                    span = n // len(doms)
                    for j, (nm, pf, hl) in enumerate(doms):
                        a, b = j * span + 3, (j + 1) * span - 3
                        dom.write(f"{sid}\t{pid}\t{n}\t{nm}\t{pf}\t{hl}\t2\t{hl - 4}\t{a}\t{b}\t"
                                  f"{r.randint(30, 300)}.5\t0.95\n")
                    if fam == "CYP":
                        cyp_tips.append((sid, pid, row["order"]))
    for f in (calls, dom, comp):
        f.close()

    # MATE below the gathering threshold: two validated copies in every 5th
    # fish, one that failed validation in every 10th. Known answers for the
    # dashboard's counting of validated below-threshold genes.
    import gzip
    os.makedirs(os.path.join(root, "rescue"), exist_ok=True)
    rc = open(os.path.join(root, "rescue_calls.tsv"), "w")
    rc.write("species\tfamily\tprotein\tprot_len\thmm_cov\tbest_evalue\tcopies\tcys_fraction\t"
             "validation\treason\n")
    added = {}
    fish = [row for row in rows if row["class"] == "Actinopteri"]
    for i, row in enumerate(fish):
        if i % 5:
            continue
        sid, tag = row["sample_id"], row["species"].split("_")[0][-6:]
        lines = []
        with open(os.path.join(root, "proteomes", sid + ".faa"), "a") as fa:
            for j, ok in enumerate([True, True] + ([False] if i % 10 == 0 else [])):
                pid = f"{tag}.mate-below{j + 1}.t1"
                fa.write(f">{pid}\n{seq_for(sid + pid, 480)}\n")
                rc.write(f"{sid}\tMATE\t{pid}\t480\t0.82\t1e-05\t{2 if ok else 1}\t\t"
                         f"{'PASS' if ok else 'FAIL'}\t{'' if ok else 'copies 1 < 2'}\n")
                if ok:
                    added[sid] = added.get(sid, 0) + 1
                    for a_, b_ in ((20, 200), (250, 440)):
                        lines.append(f"{pid} - 480 MatE PF01554.21 162 1e-05 30.1 0.1 1 2 1e-05 1e-05 "
                                     f"18.2 0.1 3 150 {a_} {b_} {a_} {b_} 0.9 -")
        with gzip.open(os.path.join(root, "rescue", sid + ".domtblout.gz"), "wt") as fh:
            fh.write("#\n" + "\n".join(lines) + "\n")
    rc.close()
    with open(os.path.join(root, "counts_rescued.tsv"), "w") as f:
        f.write("species\t" + "\t".join(FAMS) + "\n")
        for row in rows:
            v = [row["counts"][fam] + (added.get(row["sample_id"], 0) if fam == "MATE" else 0) for fam in FAMS]
            f.write(row["sample_id"] + "\t" + "\t".join(map(str, v)) + "\n")

    def nwk(tips):
        if len(tips) == 1:
            return f"{tips[0][0]}|{tips[0][1]}:0.05"
        h = len(tips) // 2
        return f"({nwk(tips[:h])},{nwk(tips[h:])})0.9:0.1"
    import json
    sys.setrecursionlimit(10000)
    with open(os.path.join(root, "trees", "CYP.tre"), "w") as f:
        f.write(nwk(cyp_tips) + ";\n")
    ins = [t for t in cyp_tips if t[2] == "Insectaformes"]
    with open(os.path.join(root, "trees", "CYP__order-insectaformes.tre"), "w") as f:
        f.write(nwk(ins) + ";\n")
    json.dump({"label": "CYP__order-insectaformes", "family": "CYP", "n_tips": len(ins),
               "n_species": len({t[0] for t in ins}), "where": {"order": ["insectaformes"]},
               "complete_only": False, "clan": [], "per_species": 0, "of_family": len(cyp_tips)},
              open(os.path.join(root, "trees", "CYP__order-insectaformes.json"), "w"))


if __name__ == "__main__":
    build(sys.argv[1], catalogue="--catalogue" in sys.argv, genes="--genes" in sys.argv)
