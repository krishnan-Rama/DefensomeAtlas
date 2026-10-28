#!/usr/bin/env bash
# Small proteins: the rescue pass, its validation rules, the metallothionein
# composition screen, and short-ORF recovery. Every number here is known.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
D="$ROOT/defensome.py"
# File modes do not always survive a tarball or a copy, and defensome.py
# executes the stub directly, so make sure it is runnable.
chmod +x "$HERE/hmmsearch_stub.py" 2>/dev/null || true

echo "== setup: extract only the defensome HMMs, and record the release =="
python3 "$HERE/make_small_protein_fixture.py" "$TMP"
export DEFENSOME_HMMSEARCH="$HERE/hmmsearch_stub.py"
python3 "$D" setup --pfam "$TMP/Pfam-A.hmm" --db-dir "$TMP/db" >/dev/null
python3 "$HERE/check_manifest.py" "$TMP/db/MANIFEST.tsv"
python3 "$D" scan --proteomes "$TMP/peps" --pfam "$TMP/db/defensome.hmm" --out "$TMP/res" >/dev/null
echo
echo "== rescue pass and validation rules =="
python3 "$D" annotate --proteomes "$TMP/peps" --out "$TMP/res" >/dev/null
python3 - "$TMP" <<'PY'
import sys, csv
t = sys.argv[1]
rc = {r["protein"]: r for r in csv.DictReader(open(f"{t}/res/rescue_calls.tsv"), delimiter="\t")}
cc = {r["protein"] for r in csv.DictReader(open(f"{t}/res/composition_candidates.tsv"), delimiter="\t")}
counts = list(csv.DictReader(open(f"{t}/res/counts.tsv"), delimiter="\t"))[0]
bad = 0
def chk(n, c):
    global bad
    print(("  OK   " if c else "  FAIL ") + n)
    if not c: bad += 1
chk("a below-threshold MATE with two repeats is rescued", rc["mate_two"]["validation"] == "PASS")
chk("  one with a single repeat fails copies>=2", rc["mate_one"]["validation"] == "FAIL")
chk("a 43 aa Cys-rich hit is rescued as metallothionein", rc["mt_real_1"]["validation"] == "PASS")
chk("  a 300 aa spurious hit fails len<=120", rc["mt_fp_long"]["validation"] == "FAIL")
chk("rescued calls stay OUT of counts.tsv", counts["MATE"] == "0" and counts["Metallothionein"] == "0")
chk("composition screen keeps both real MTs", {"mt_real_1", "mt_real_2"} <= cc)
chk("  and rejects the defensin decoy (6 Cys, aromatic)", "defensin" not in cc)
chk("  and rejects a Cys-rich but aromatic peptide", "cys_arom" not in cc)
if bad: sys.exit(1)
PY

echo
echo "== qc verdicts: a zero is not automatically a bad map rule =="
python3 "$D" qc --out "$TMP/res" > "$TMP/qc.log" 2>&1
python3 - "$TMP" <<'PY'
import sys, csv
t = sys.argv[1]
z = {r["family"]: r for r in csv.DictReader(open(f"{t}/res/qc_zero_families.tsv"), delimiter="\t")}
f = {r[""]: r for r in csv.DictReader(open(f"{t}/res/qc_families.tsv"), delimiter="\t")}
bad = 0
def chk(n, c):
    global bad
    print(("  OK   " if c else "  FAIL ") + n)
    if not c: bad += 1
chk("MATE reported as found below the gathering threshold", z["MATE"]["verdict"] == "FOUND_BELOW_GA")
chk("metallothionein likewise, despite a sibling accession being absent",
    z["Metallothionein"]["verdict"] == "FOUND_BELOW_GA"
    and "PF00131" in z["Metallothionein"]["accessions_missing"]
    and "only partly searched" in z["Metallothionein"]["what_it_means"])
chk("a family with NO accession in the database says so",
    any(v["verdict"] == "ACCESSION_NOT_IN_DATABASE" for v in z.values()))
chk("control families are NOT blamed on the map rule when never searched",
    all(v["verdict"] != "MAP_RULE_TOO_STRICT" for v in f.values()))
if bad: sys.exit(1)
PY

echo
echo "== short-ORF recovery from transcripts =="
python3 "$HERE/make_orf_fixture.py" "$TMP/orf"
python3 "$D" short-orfs --samplesheet "$TMP/orf/samplesheet.tsv" --out "$TMP/orfres" >/dev/null
python3 - "$TMP" <<'PY'
import sys, csv
t = sys.argv[1]
h = list(csv.DictReader(open(f"{t}/orfres/short_orf_mt.tsv"), delimiter="\t"))
s = list(csv.DictReader(open(f"{t}/orfres/short_orf_summary.tsv"), delimiter="\t"))
g = {r["gene"] for r in h}
bad = 0
def chk(n, c):
    global bad
    print(("  OK   " if c else "  FAIL ") + n)
    if not c: bad += 1
chk("metallothionein found on the forward strand",
    any(r["gene"] == "TRINITY_DN10_c0_g1" and r["strand"] == "+" for r in h))
chk("and one on the reverse strand",
    any(r["gene"] == "TRINITY_DN22_c0_g1" and r["strand"] == "-" for r in h))
chk("the recovered sequence is exactly Drosophila MtnB",
    "MVCKGCGTNCQCSAQKCGDNCACNKDCQCVCKNGPKDQCCSNK" in {r["protein"] for r in h})
chk("two isoforms of one gene count as one gene", s[0]["mt_like_genes"] == "2")
chk("the defensin decoy is rejected", "TRINITY_DN33_c0_g1" not in g)
chk("an ORF above the length band is excluded", "TRINITY_DN44_c0_g1" not in g)
chk("a sample with no transcripts is skipped, not failed", len(s) == 1)
if bad: sys.exit(1)
PY
echo
echo "all small-protein tests passed"
