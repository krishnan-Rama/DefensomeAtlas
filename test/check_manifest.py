#!/usr/bin/env python3
"""The manifest is what lets qc say an accession was never searched."""
import csv, sys
m = list(csv.DictReader(open(sys.argv[1]), delimiter="\t"))
found = [r for r in m if r["found"] == "yes"]
assert len(m) > 30, f"manifest has only {len(m)} rows"
assert found, "no accession found in the database"
print(f"  OK   manifest lists all {len(m)} map accessions; {len(found)} present in this release")
