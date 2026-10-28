#!/usr/bin/env python3
"""Add a `pipeline` column to an existing sample sheet from its FASTA IDs."""
import re, sys
import pandas as pd
HINTS = [("BRAKER", r"^BRAKER"), ("Ensembl", r"^ENS[A-Z]*P\d"),
         ("AUGUSTUS", r"^AUGUSTUS"), ("Trinity", r"_c\d+_g\d+_i\d+")]
p = sys.argv[1]
t = pd.read_csv(p, sep="\t")
def first_ids(fa, k=200):
    out = []
    with open(fa) as fh:
        for l in fh:
            if l.startswith(">"):
                out.append(l[1:].split()[0])
                if len(out) >= k: break
    return out
pipes = []
for fa in t.fasta:
    ids = first_ids(fa)
    pipes.append(next((n for n, rx in HINTS
                       if sum(1 for i in ids if re.search(rx, i)) > len(ids) * .5), "unknown"))
t.insert(3, "pipeline", pipes)
t.to_csv(p, sep="\t", index=False)
print(t[["sample_id", "method", "pipeline", "est_genes"]].to_string(index=False))
