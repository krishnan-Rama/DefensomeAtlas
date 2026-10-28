#!/usr/bin/env python3
"""Build a tiny multi-method proteome set with a known answer.

3 species x 3 methods. Every transcriptome sample has exactly 120 genes spread
over a variable number of isoforms; every DToL sample has exactly 150 genes and
no isoforms. CYP counts are injected at 20 / 16 / 12 so recovery is exactly
0.80 and 0.60. GST is injected as 8 real N+C pairs plus 6 unpaired GST_N
decoys, which the ALL rule must reject. utrorf entries are only ever placed on
a non-first isoform, so dropping them removes a protein but never a gene."""
import gzip, os, random, sys

SPECIES = ["X_alpha", "Y_beta", "Z_gamma"]
CYP_N = {"DToL": 20, "genome_guided": 16, "denovo": 12}


def aa(n):
    return "M" + "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n - 1))


def build(root):
    random.seed(1)
    p = os.path.join(root, "peps"); os.makedirs(p, exist_ok=True)
    for sp in SPECIES:
        with open(f"{p}/{sp}_GG_010101_okay.aa.fasta", "w") as f:
            for g in range(120):
                for i in range(1, random.choice([1, 1, 2, 3]) + 1):
                    f.write(f">{sp}_010101_GG_{g//20}_c{g%20}_g1_i{i}\n{aa(random.randint(150,600))}\n")
        with open(f"{p}/{sp}_020202_okay.aa.fasta", "w") as f:
            for g in range(120):
                n_iso = random.choice([2, 3, 6])   # always >1 so a dropped
                for i in range(1, n_iso + 1):      # utrorf cannot delete a gene
                    u = "utrorf" if (g % 17 == 0 and i == n_iso) else ""
                    f.write(f">{sp}_020202_{g//20}_c{g%20}_g1_i{i}{u}\n{aa(random.randint(150,600))}\n")
        with open(f"{p}/{sp}_pep_DToL.fasta", "w") as f:
            for g in range(150):
                f.write(f">{sp}P{g:06d}.1 pep gene:{sp}G{g:06d}\n{aa(random.randint(200,600))}\n")


def domtbl(root):
    random.seed(2)
    res = os.path.join(root, "res")
    hd = os.path.join(res, "hmmsearch"); os.makedirs(hd, exist_ok=True)
    for fn in sorted(os.listdir(os.path.join(res, "proteomes"))):
        sid = fn[:-4]
        meth = "DToL" if "DToL" in sid else ("genome_guided" if "_GG_" in sid else "denovo")
        ids = [l[1:].split()[0] for l in open(os.path.join(res, "proteomes", fn))
               if l.startswith(">")]
        out, k = [], 0

        def row(pid, tl, nm, pf, hl, cov=0.9):
            return (f"{pid} - {tl} {nm} {pf}.1 {hl} 1e-50 200 0.1 1 1 1e-50 1e-50 "
                    f"200 0.1 1 {int(hl*cov)} 1 {tl} 1 {tl} 0.95 d\n")

        for _ in range(CYP_N[meth]):
            out.append(row(ids[k], 520, "p450", "PF00067", 480)); k += 1
        for _ in range(8):                       # complete GSTs
            out.append(row(ids[k], 230, "GST_N", "PF02798", 80))
            out.append(row(ids[k], 230, "GST_C", "PF00043", 100)); k += 1
        for _ in range(6):                       # unpaired decoys, must be rejected
            out.append(row(ids[k], 240, "GST_N", "PF02798", 80)); k += 1
        out.append(row(ids[k], 700, "Glu_cys_ligase", "PF03074", 640)); k += 1
        with gzip.open(f"{hd}/{sid}.domtblout.gz", "wt") as f:
            f.write("# fixture\n"); f.writelines(out)


if __name__ == "__main__":
    root = sys.argv[1]
    if "--domtbl" in sys.argv:
        domtbl(root)
    else:
        build(root)
