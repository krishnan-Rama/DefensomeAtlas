#!/usr/bin/env bash
# The dashboard as a site, and subset trees. Needs no MAFFT or FastTree: stubs
# stand in for them, because what is tested is which sequences go in and what
# the outputs are called, not the phylogenetics.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
D="python3 $ROOT/defensome.py"
python3 "$HERE/make_taxonomy_fixture.py" "$TMP/res" --genes >/dev/null
R=$TMP/res

echo "== extract, then trees on subsets (stub aligner and tree builder) =="
$D extract --out "$R" --proteomes "$R/proteomes" >"$TMP/x.log" 2>&1 || { cat "$TMP/x.log"; exit 1; }
mkdir -p "$TMP/bin"
cat > "$TMP/bin/mafft" <<'SH'
#!/usr/bin/env bash
cat "${@: -1}"
SH
cat > "$TMP/bin/FastTree" <<'SH'
#!/usr/bin/env bash
python3 - "${@: -1}" <<'PY'
import sys
ids = [l[1:].split()[0] for l in open(sys.argv[1]) if l.startswith(">")]
print("(" + ",".join(f"{i}:0.1" for i in ids) + ");")
PY
SH
chmod +x "$TMP/bin/"*
export PATH="$TMP/bin:$PATH"
$D trees --out "$R" --families GST --metadata "$R/metadata.tsv" --where class=Insecta \
   >"$TMP/t1.log" 2>&1 || { cat "$TMP/t1.log"; exit 1; }
$D trees --out "$R" --families GST --metadata "$R/metadata.tsv" --where "habitat=marine,freshwater" \
   --where phylum=Chordata --complete-only >"$TMP/t2.log" 2>&1 || { cat "$TMP/t2.log"; exit 1; }
if $D trees --out "$R" --families GST --metadata "$R/metadata.tsv" --where diet=x >"$TMP/t3.log" 2>&1; then
  echo "  FAIL an unknown --where column was accepted"; exit 1; fi
python3 - "$R" <<'PY'
import sys, json, csv, os
R = sys.argv[1]
md = {r["sample_id"]: r for r in csv.DictReader(open(f"{R}/metadata.tsv"), delimiter="\t")}
calls = list(csv.DictReader(open(f"{R}/gene_calls_all.tsv"), delimiter="\t"))
comp = {(r["species"], r["protein"]): r["status"]
        for r in csv.DictReader(open(f"{R}/domains/completeness.tsv"), delimiter="\t") if r["family"] == "GST"}
def tips(label):
    t = open(f"{R}/trees/{label}.tre").read().strip().strip("();")
    return {x.rsplit(":", 1)[0] for x in t.split(",")}
bad = 0
def chk(n, c):
    global bad
    print(("  OK   " if c else "  FAIL ") + n); bad += not c
want1 = {f"{r['species']}|{r['protein']}" for r in calls
         if r["family"] == "GST" and md[r["species"]]["class"] == "Insecta"}
got1 = tips("GST__class-insecta")
chk(f"--where class=Insecta keeps exactly the insect GSTs ({len(got1)})", got1 == want1)
m1 = json.load(open(f"{R}/trees/GST__class-insecta.json"))
chk("the sidecar records the restriction and the counts",
    m1["where"] == {"class": ["insecta"]} and m1["n_tips"] == len(want1) and m1["n_species"] == 46
    and m1["of_family"] == sum(r["family"] == "GST" for r in calls))
want2 = {f"{r['species']}|{r['protein']}" for r in calls if r["family"] == "GST"
         and md[r["species"]]["habitat"] in ("marine", "freshwater")
         and md[r["species"]]["phylum"] == "Chordata" and comp[(r["species"], r["protein"])] == "COMPLETE"}
lab2 = "GST__habitat-marine+freshwater__phylum-chordata__complete"
got2 = tips(lab2)
chk(f"two --where (AND, with OR inside) plus --complete-only ({len(got2)})", got2 == want2 and len(want2) > 0)
chk("the whole-family tree name is never reused for a subset", not os.path.exists(f"{R}/trees/GST.tre"))
sys.exit(1 if bad else 0)
PY
grep -q "is not in" "$TMP/t3.log" && echo "  OK   an unknown --where column stops with the column list"

echo "== --thin, the length filter, and a MAFFT that runs out of memory =="
python3 - "$R" <<'PY'
import sys
R = sys.argv[1]
# a read-through model in one insect: four GSTs end to end
with open(f"{R}/fasta/GST.faa", "a") as fh:
    fh.write(">Genusins070_species070_pep|ins070.readthrough.t1\n" + "MKLV" * 500 + "\n")
PY
$D trees --out "$R" --families GST --metadata "$R/metadata.tsv" --where phylum=Arthropoda \
   --thin class:2 >"$TMP/t4.log" 2>&1 || { cat "$TMP/t4.log"; exit 1; }
$D trees --out "$R" --families GST --metadata "$R/metadata.tsv" --where class=Insecta --force \
   >"$TMP/t5.log" 2>&1 || { cat "$TMP/t5.log"; exit 1; }
python3 - "$R" <<'PY'
import sys, json, csv
R = sys.argv[1]
md = {r["sample_id"]: r for r in csv.DictReader(open(f"{R}/metadata.tsv"), delimiter="\t")}
calls = [r for r in csv.DictReader(open(f"{R}/gene_calls_all.tsv"), delimiter="\t") if r["family"] == "GST"]
def tips(label):
    t = open(f"{R}/trees/{label}.tre").read().strip().strip("();")
    return {x.rsplit(":", 1)[0] for x in t.split(",")}
bad = 0
def chk(n, c):
    global bad
    print(("  OK   " if c else "  FAIL ") + n); bad += not c
got = tips("GST__phylum-arthropoda__thin-class2")
sps = {t.split("|")[0] for t in got}
by = {}
for sp in sps: by.setdefault(md[sp]["class"], []).append(sp)
chk(f"--thin class:2 keeps two species per class ({sorted((k, len(v)) for k, v in by.items())})",
    sorted(by) == ["Insecta", "Malacostraca"] and all(len(v) == 2 for v in by.values()))
n = {}
for r in calls: n[r["species"]] = n.get(r["species"], 0) + 1
chk("each kept species keeps every one of its copies",
    all(sum(t.startswith(sp + "|") for t in got) == n[sp] for sp in sps))
ins = [sp for sp in md if md[sp]["class"] == "Insecta"]
med = sorted(n[s] for s in ins)[len(ins) // 2]
chk("and they are the species with the most typical copy number",
    all(abs(n[s] - med) <= min(abs(n[x] - med) for x in ins if x not in sps) for s in by["Insecta"]))
chk("a read-through model is dropped by the length filter",
    not any("readthrough" in t for t in tips("GST__class-insecta")))
m = json.load(open(f"{R}/trees/GST__class-insecta.json"))
chk("and the sidecar says so", m["dropped_long"] == 1)
sys.exit(1 if bad else 0)
PY
grep -q "dropped (longest 2,000 aa)" "$TMP/t5.log" && echo "  OK   the log names the dropped length"
cat > "$TMP/bin/mafft" <<'SH'
#!/usr/bin/env bash
echo "line 2842: 1012019 Killed  disttbfast" >&2; exit 1
SH
if $D trees --out "$R" --families GST,CYP --metadata "$R/metadata.tsv" --where class=Mammalia \
     >"$TMP/t6.log" 2>&1 && grep -q "out of memory" "$TMP/t6.log" && ! grep -q Traceback "$TMP/t6.log" \
     && [[ $(grep -c "MAFFT failed" "$TMP/t6.log") == 2 ]]; then
  echo "  OK   a killed MAFFT is reported as out of memory and the next family still runs"
else echo "  FAIL"; cat "$TMP/t6.log"; exit 1; fi
cat > "$TMP/bin/mafft" <<'SH'
#!/usr/bin/env bash
cat "${@: -1}"
SH

echo "== dashboard: sequences from fasta/ when --proteomes is not given =="
$D dashboard --out "$R" --metadata "$R/metadata.tsv" --group-by class >"$TMP/d0.log" 2>&1 || { cat "$TMP/d0.log"; exit 1; }
grep -q "sequences from" "$TMP/d0.log" && echo "  OK   falls back to the extracted family FASTA"
grep -q "0 not found" "$TMP/d0.log" && echo "  OK   every called protein found" || { echo "  FAIL"; grep sequences "$TMP/d0.log"; exit 1; }

echo "== structure: Boltz input for a CYP with its haem and a ligand =="
$D structure --out "$R" --gene 'Genusins070_species070_pep|ins070.cyp-g3.t1' \
   --ligand 'CC1(C)C(C=C(Cl)Cl)C1C(=O)OCc1cccc(Oc2ccccc2)c1' --ligand-name permethrin \
   --proteomes "$R/proteomes" >"$TMP/s1.log" 2>&1 || { cat "$TMP/s1.log"; exit 1; }
python3 - "$R" <<'PY'
import sys, glob, json
R = sys.argv[1]
d = glob.glob(f"{R}/structures/*permethrin")[0]
y = open(glob.glob(d + "/*.yaml")[0]).read()
seq = {}
for l in open(f"{R}/proteomes/Genusins070_species070_pep.faa"):
    if l.startswith(">"): k = l[1:].split()[0]
    else: seq[k] = l.strip()
ok = ("ccd: HEM" in y and "binder: C" in y and seq["ins070.cyp-g3.t1"] in y
      and "boltz predict" in open(d + "/run_boltz.sh").read() and "--no_kernels" in open(d + "/run_boltz.sh").read())
print(("  OK   " if ok else "  FAIL ") + "the YAML holds the sequence, haem, the ligand and its affinity request")
sys.exit(0 if ok else 1)
PY
if $D structure --out "$R" --gene 'Genusins070_species070_pep|nope' >"$TMP/s2.log" 2>&1; then
  echo "  FAIL a gene that is not a call was accepted"; exit 1; fi
echo "  OK   a protein that is not a defensome call is refused"
python3 "$HERE/make_fake_boltz.py" "$(ls -d "$R"/structures/*permethrin)" >/dev/null

echo "== dashboard with --proteomes, trees loaded on demand =="
DEFENSOME_TREE_EMBED_MAX=20000 $D dashboard --out "$R" --metadata "$R/metadata.tsv" --group-by class \
  --proteomes "$R/proteomes" >"$TMP/d1.log" 2>&1 || { cat "$TMP/d1.log"; exit 1; }
ls "$R/dashboard/data/trees/" | grep -q "GST__class-insecta.js" && echo "  OK   large trees are written to data/trees"
DEFENSOME_TREE_EMBED_MAX=20000 $D dashboard --out "$R" --metadata "$R/metadata.tsv" --group-by class \
  --proteomes "$R/proteomes" >"$TMP/d2.log" 2>&1
grep -q "unchanged since the last export" "$TMP/d2.log" && echo "  OK   an unchanged rerun reuses the species data"
grep -q "structures: 1 predicted" "$TMP/d2.log" && echo "  OK   the finished prediction is picked up"
$D dashboard --out "$R" --dashboard-out "$TMP/ga" --ga-only --metadata "$R/metadata.tsv" >"$TMP/ga.log" 2>&1 || { cat "$TMP/ga.log"; exit 1; }
python3 - "$TMP/ga/index.html" "$R/counts.tsv" <<'PY'
import sys, re, json, csv
P = json.loads(re.search(r'const D = (\{.*?\});\n', open(sys.argv[1]).read(), re.S).group(1))
c = {r["species"]: int(r["MATE"]) for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
ok = [P["counts"]["MATE"][i] == c[s["name"]] for i, s in enumerate(P["species"])]
print(("  OK   " if all(ok) else "  FAIL ") + "--ga-only counts exactly counts.tsv")
sys.exit(0 if all(ok) else 1)
PY
bash "$HERE/test_dashboard.sh" "$R/dashboard/index.html" > "$TMP/smoke.log" 2>&1 \
  || { echo "  FAIL headless smoke test"; grep FAIL "$TMP/smoke.log"; exit 1; }
echo "  OK   headless smoke test"
if node -e "require(process.env.PLAYWRIGHT_MODULE || 'playwright')" 2>/dev/null; then
  echo; echo "== real browser =="
  node "$HERE/site_browser_test.js" "$R/dashboard" "$R"
else
  echo "playwright not installed; skipping the real-browser site tests"
fi

echo; echo "== --no-sequences =="
$D dashboard --out "$R" --dashboard-out "$TMP/ns" --no-sequences >"$TMP/d3.log" 2>&1 || { cat "$TMP/d3.log"; exit 1; }
python3 - "$TMP/ns" <<'PY'
import sys, glob, json
f = sorted(glob.glob(sys.argv[1] + "/data/species/*.js"))[0]
s = open(f).read(); o = json.loads(s[s.index("{"):s.rindex(")")])
ok = o["genes"] and o["domains"] and not o["seq"] and not glob.glob(sys.argv[1] + "/data/fasta/*.gz")
print(("  OK   " if ok else "  FAIL ") + "genes and domains kept, sequences and FASTA left out")
sys.exit(0 if ok else 1)
PY
echo; echo "all site tests passed"
