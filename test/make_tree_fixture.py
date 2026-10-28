#!/usr/bin/env python3
"""Known-answer fixture for the dashboard tree tab.

3 species x 3 methods. The CYP gene tree has 12 orthologous gene clades of 9
tips (every species, every method) at support 0.95, plus three DENOVO-ONLY
clades of 5, 4 and 3 tips at support 0.61 whose proteins are marked FRAGMENT.
The single-source finder must report exactly those three and nothing else."""
import os, random, sys

def build(root):
    random.seed(4)
    for d in ("trees", "cyp", "domains"):
        os.makedirs(os.path.join(root, d), exist_ok=True)
    species = ["X_alpha", "Y_beta", "Z_gamma"]
    meth = {"DToL": "pep_DToL", "genome_guided": "GG_010101_okay", "denovo": "020202_okay"}
    samples = {(sp, m): f"{sp}_{suf}" for sp in species for m, suf in meth.items()}
    with open(f"{root}/samplesheet.tsv", "w") as f:
        f.write("sample_id\tspecies\tmethod\tpipeline\tfasta\tgene_regex\n")
        for (sp, m), sid in samples.items():
            f.write(f"{sid}\t{sp}\t{m}\t{'Ensembl' if m == 'DToL' else 'Trinity'}\t/x/{sid}.faa\t\n")
    fams = ["CYP", "CCE", "GST", "UGT", "ABC_full", "GCL", "Catalase"]
    with open(f"{root}/counts.tsv", "w") as f:
        f.write("species\t" + "\t".join(fams) + "\n")
        for sid in samples.values():
            f.write(sid + "\t" + "\t".join("1" if x == "GCL" else str(random.randint(5, 60)) for x in fams) + "\n")
    clans = ["CYP2", "CYP3", "CYP4", "MITO"]
    calls, genes = [], []
    for g in range(12):
        tips = []
        for sp in species:
            for m in meth:
                sid, pid = samples[(sp, m)], f"{sp}_{m}_g{g}"
                tips.append(f"{sid}|{pid}:0.0{random.randint(1, 9)}")
                calls.append((sid, pid, clans[g % 4]))
        genes.append("(" + ",".join(tips) + ")0.95:0.2")
    for j, n in enumerate([4, 3, 5]):
        tips = []
        for k in range(n):
            sp = species[k % 3]; sid = samples[(sp, "denovo")]; pid = f"{sp}_denovo_frag{j}_{k}"
            tips.append(f"{sid}|{pid}:0.03"); calls.append((sid, pid, "CYP3"))
        genes.append("(" + ",".join(tips) + ")0.61:0.3")
    open(f"{root}/trees/CYP.tre", "w").write("(" + ",".join(genes) + ");\n")
    with open(f"{root}/cyp/CYP_clan_calls.tsv", "w") as f:
        f.write("species\tprotein\tclan\tscore\tmargin\tconfidence\n")
        for sid, pid, c in calls:
            f.write(f"{sid}\t{pid}\t{c}\t500\t300\tOK\n")
    with open(f"{root}/domains/completeness.tsv", "w") as f:
        f.write("species\tprotein\tfamily\tcategory\tstatus\tn_required\tn_complete\t"
                "min_domain_cov\tmean_domain_cov\tprot_len\tlength_ratio\tdomain_coverages\n")
        for sid, pid, c in calls:
            st = "FRAGMENT" if "frag" in pid else "COMPLETE"
            f.write(f"{sid}\t{pid}\tCYP\tPhaseI\t{st}\t1\t1\t0.9\t0.9\t500\t1.0\tPF00067:0.90\n")

if __name__ == "__main__":
    build(sys.argv[1])
