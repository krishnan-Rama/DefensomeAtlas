#!/usr/bin/env python3
"""Write a metadata template from a sample sheet, for between-species grouping.

`compare` pairs annotations of one species. When you have one proteome per
species, the contrast is between species instead, and that needs a grouping
variable you supply: lineage, habitat, exposure, whatever the question is.
This writes the skeleton with the species already filled in."""
import csv, sys, os

if len(sys.argv) < 3:
    sys.exit("usage: make_metadata_template.py samplesheet.tsv metadata.tsv [col ...]")
sheet, out = sys.argv[1], sys.argv[2]
cols = sys.argv[3:] or ["lineage", "habitat", "notes"]
if os.path.exists(out):
    sys.exit(f"ERROR: {out} exists; refusing to overwrite it")
rows = list(csv.DictReader(open(sheet), delimiter="\t"))
seen, order = set(), []
for r in rows:
    s = (r.get("species") or r["sample_id"]).strip()
    if s not in seen:
        seen.add(s); order.append((r["sample_id"], s, (r.get("pipeline") or "").strip()))
with open(out, "w") as f:
    f.write("\t".join(["sample_id", "species", "pipeline"] + cols) + "\n")
    for sid, sp, pipe in order:
        f.write("\t".join([sid, sp, pipe] + [""] * len(cols)) + "\n")
print(f"wrote {out} with {len(order)} rows and blank columns: {', '.join(cols)}")
print("Fill them in, then:")
print(f"  python3 defensome.py dashboard --out <results> --metadata {out} --group-by {cols[0]}")
print("The grouping is yours to choose and to justify; nothing here can guess it.")
