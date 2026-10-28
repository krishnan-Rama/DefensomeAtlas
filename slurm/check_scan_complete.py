#!/usr/bin/env python3
"""List sample_ids from a sample sheet that have no usable scan output.

    check_scan_complete.py samplesheet.tsv results_dir

Prints one missing sample_id per line; prints nothing when the scan is
complete. A zero-byte file counts as missing: an interrupted task can leave
one behind, and annotating it would silently drop that species."""
import csv, os, sys

if len(sys.argv) < 3:
    sys.exit(__doc__)
sheet, out = sys.argv[1], sys.argv[2]
miss = []
for r in csv.DictReader(open(sheet), delimiter="\t"):
    sid = r["sample_id"]
    p = os.path.join(out, "hmmsearch", sid + ".domtblout.gz")
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        miss.append(sid)
print("\n".join(miss))
