#!/usr/bin/env bash
# End-to-end test that needs no external tools. Builds a small multi-method
# proteome set with a KNOWN answer, runs inspect -> collapse -> annotate ->
# compare, and checks the numbers come back right.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
python3 "$HERE/make_fixture.py" "$TMP"

echo "== inspect =="
python3 "$ROOT/defensome.py" inspect --proteomes "$TMP/peps" --out "$TMP" >"$TMP/inspect.log"
grep -q "trinity_isoform" "$TMP/inspect.log" || { echo "FAIL: Trinity rule not detected"; exit 1; }
grep -q "gene_tag"        "$TMP/inspect.log" || { echo "FAIL: gene: tag not detected"; exit 1; }
echo "  OK   gene-ID rules detected"

echo "== collapse =="
python3 "$ROOT/defensome.py" collapse --samplesheet "$TMP/samplesheet.tsv" --out "$TMP/res" \
        --drop-utrorf >"$TMP/collapse.log"
python3 - "$TMP" <<'PY'
import sys, csv
t=sys.argv[1]
rows=list(csv.DictReader(open(f"{t}/res/collapse_stats.tsv"), delimiter="\t"))
by={}
for r in rows: by.setdefault(r["method"],[]).append(r)
# the fixture builds exactly 120 genes per transcriptome sample and 150 per DToL
for m,rs in by.items():
    for r in rs:
        want = 150 if m=="DToL" else 120
        got = int(r["n_genes_out"])
        assert got==want, f"FAIL {r['sample_id']}: {got} genes, expected {want}"
assert any(int(r["dropped_utrorf"])>0 for r in rows), "FAIL: no utrorf dropped"
print("  OK   isoform collapse exact for every sample")
print("  OK   utrorf entries dropped")
PY

echo "== annotate =="
python3 "$HERE/make_fixture.py" "$TMP" --domtbl
python3 "$ROOT/defensome.py" annotate --proteomes "$TMP/res/proteomes" --out "$TMP/res" \
        >"$TMP/annotate.log" 2>&1
python3 - "$TMP" <<'PY'
import sys, csv
t=sys.argv[1]
r=list(csv.DictReader(open(f"{t}/res/counts.tsv"), delimiter="\t"))
assert len(r)==9, f"FAIL: {len(r)} samples, expected 9"
# the fixture injects 20 CYPs in DToL, 16 in GG, 12 in denovo
for row in r:
    sid=row[list(row)[0]]
    want = 20 if "DToL" in sid else (16 if "_GG_" in sid else 12)
    assert int(row["CYP"])==want, f"FAIL {sid}: CYP={row['CYP']}, expected {want}"
    assert int(row["GST"])==8, f"FAIL {sid}: GST={row['GST']}, expected 8 (ALL rule)"
print("  OK   per-method CYP counts exact")
print("  OK   GST ALL-rule rejects unpaired GST_N decoys")
PY

echo "== compare =="
python3 "$ROOT/defensome.py" compare --samplesheet "$TMP/samplesheet.tsv" \
        --out "$TMP/res" --reference DToL >"$TMP/compare.log" 2>&1
python3 - "$TMP" <<'PY'
import sys, csv
t=sys.argv[1]
rec={}
for r in csv.DictReader(open(f"{t}/res/comparison/recovery.tsv"), delimiter="\t"):
    if r["family"]=="CYP": rec.setdefault(r["method"],set()).add(round(float(r["recovery"]),2))
assert rec["genome_guided"]=={0.8}, f"FAIL: GG CYP recovery {rec['genome_guided']}, expected 0.80"
assert rec["denovo"]=={0.6},        f"FAIL: denovo CYP recovery {rec['denovo']}, expected 0.60"
print("  OK   recovery matches the injected truth exactly (GG 0.80, denovo 0.60)")
PY
echo "== compare --exclude =="
python3 "$ROOT/defensome.py" compare --samplesheet "$TMP/samplesheet.tsv" \
        --out "$TMP/res" --exclude X_alpha >"$TMP/compare_ex.log" 2>&1
grep -q "2 of 2 species have every method" "$TMP/compare_ex.log" \
  || { echo "FAIL: --exclude did not drop a species"; exit 1; }
echo "  OK   --exclude drops a species from the paired design"
echo; echo "all pipeline tests passed"
