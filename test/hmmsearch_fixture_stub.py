#!/usr/bin/env python3
"""hmmsearch stand-in for the gene-level fixture. A protein named
<tag>.<family>-g<n>.t1 hits every Pfam model of that family over its whole
length, at the gathering threshold only; nothing else hits. So the right
answer for any proteome is computable from its headers and the map."""
import csv, os, sys
a = sys.argv
out = a[a.index("--domtblout") + 1]; fa = a[-1]
here = os.path.dirname(os.path.abspath(__file__))
mp = os.environ.get("DEFENSOME_TEST_MAP", os.path.join(here, "..", "defensome_map.tsv"))
fam = {}
for r in csv.DictReader((l for l in open(mp) if not l.startswith("#")), delimiter="\t"):
    fam[r["family"].lower()] = [x.strip() for x in r["pfam_ids"].split(",")]
L, cur = {}, None
for line in open(fa):
    if line.startswith(">"):
        cur = line[1:].split()[0]; L[cur] = 0
    elif cur:
        L[cur] += len(line.strip())
rows = []
if "--cut_ga" in a:
    for pid, n in L.items():
        if "-g" not in pid or "." not in pid:
            continue
        f = pid.split(".", 1)[1].rsplit("-g", 1)[0]
        for k, acc in enumerate(fam.get(f, [])):
            a_, b_ = (1, n) if len(fam[f]) == 1 else ((1, n // 2) if k == 0 else (n // 2 + 1, n))
            rows.append(f"{pid} - {n} M{acc} {acc}.1 100 1e-30 200 0.1 1 1 1e-30 1e-30 200 0.1 "
                        f"1 100 {a_} {b_} {a_} {b_} 0.99 stub\n")
with open(out, "w") as fh:
    fh.write("# stub\n"); fh.writelines(rows)
