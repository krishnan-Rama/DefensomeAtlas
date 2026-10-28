#!/usr/bin/env python3
# GA pass: only the CYP hits. Rescue pass (no --cut_ga): sub-threshold hits.
import sys, os
a = sys.argv
out = a[a.index("--domtblout") + 1]; fa = a[-1]
L = {}; cur = None
for line in open(fa):
    if line.startswith(">"): cur = line[1:].split()[0]; L[cur] = 0
    else: L[cur] += len(line.strip())
def row(t, qname, qacc, qlen, hf, ht, af, at, sc):
    return (f"{t} - {L[t]} {qname} {qacc}.1 {qlen} 1e-3 {sc} 0.1 1 1 1e-3 1e-3 {sc} 0.1 "
            f"{hf} {ht} {af} {at} {af} {at} 0.9 d\n")
rows = []
if "--cut_ga" in a:
    if "cyp_ga" in L: rows.append(row("cyp_ga","p450","PF00067",459,1,440,20,480,300))
else:
    if "mt_real_1" in L:  rows.append(row("mt_real_1","Metallothio_5","PF02067",42,1,40,2,42,12))
    if "mt_fp_long" in L: rows.append(row("mt_fp_long","Metallothio_5","PF02067",42,2,38,100,140,11))
    if "mate_two" in L:
        rows.append(row("mate_two","MatE","PF01554",161,1,158,40,200,18))
        rows.append(row("mate_two","MatE","PF01554",161,2,159,290,450,17))
    if "mate_one" in L:   rows.append(row("mate_one","MatE","PF01554",161,1,157,40,200,16))
with open(out, "w") as f:
    f.write("# stub\n"); f.writelines(rows)
