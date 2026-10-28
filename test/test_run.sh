#!/usr/bin/env bash
# `run` end to end on a chosen subset of species, with the prebuilt HMM
# database fetched from a release-style bundle. hmmsearch is a stub whose
# answer is computable from the headers, so the counts are checked exactly.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
D="python3 $ROOT/defensome.py"
python3 "$HERE/make_taxonomy_fixture.py" "$TMP/fx" --genes >/dev/null
export DEFENSOME_HMMSEARCH="$HERE/hmmsearch_fixture_stub.py"
chmod +x "$DEFENSOME_HMMSEARCH"

echo "== a release bundle built by setup --bundle, fetched by run =="
python3 - "$ROOT/defensome_map.tsv" "$TMP/Pfam-A.hmm" <<'PY'
import sys, csv
accs = sorted({x.strip() for r in csv.DictReader((l for l in open(sys.argv[1]) if not l.startswith("#")),
               delimiter="\t") for x in r["pfam_ids"].split(",")})
with open(sys.argv[2], "w") as f:
    for a in accs + ["PF99999"]:            # one model the map does not use
        f.write(f"HMMER3/f [3.4]\nNAME  M{a}\nACC   {a}.1\nDESC  stub\nLENG  100\nGA    20.0 20.0;\nHMM\n//\n")
PY
$D setup --pfam "$TMP/Pfam-A.hmm" --db-dir "$TMP/rel" --pfam-version "Pfam release 99.0 (stub)" \
   --bundle >"$TMP/setup.log" 2>&1 || { cat "$TMP/setup.log"; exit 1; }
[[ -s "$TMP/rel/defensome_db.tar.gz" && -s "$TMP/rel/defensome_db.tar.gz.sha256" ]] \
  && echo "  OK   setup --bundle writes the release asset and its checksum"
if grep -q "PF99999" <(gzip -dc "$TMP/rel/defensome_db.tar.gz" | strings); then
  echo "  FAIL the bundle carries a model the map does not use"; exit 1; fi
echo "  OK   the bundle holds only the map's models"

printf 'Genusins070 species070\nGenusact001_species001_pep\n# a comment\nGenusbiv121\n' > "$TMP/species.txt"
$D run --proteomes "$TMP/fx/proteomes" --out "$TMP/out" --species "$TMP/species.txt" \
   --metadata "$TMP/fx/metadata.tsv" --group-by class --db-url "$TMP/rel/defensome_db.tar.gz" \
   --threads 2 >"$TMP/run.log" 2>&1 || { tail -40 "$TMP/run.log"; exit 1; }
python3 - "$TMP" "$ROOT/defensome_map.tsv" <<'PY'
import sys, csv, os, glob
T, mp = sys.argv[1], sys.argv[2]
rules = {r["family"]: (int(r["min_len"]), r["rule"], {x.strip() for x in r["pfam_ids"].split(",")})
         for r in csv.DictReader((l for l in open(mp) if not l.startswith("#")), delimiter="\t")}
bad = 0
def chk(n, c):
    global bad
    print(("  OK   " if c else "  FAIL ") + n); bad += not c
counts = {r["species"]: r for r in csv.DictReader(open(f"{T}/out/counts.tsv"), delimiter="\t")}
want_sp = {"Genusins070_species070_pep", "Genusact001_species001_pep", "Genusbiv121_species121_pep"}
chk(f"three names (binomial, file stem, genus) select three proteomes ({sorted(counts)})", set(counts) == want_sp)
ok = True
for sp in want_sp:
    exp = {}
    L, cur = {}, None
    for line in open(f"{T}/fx/proteomes/{sp}.faa"):
        if line.startswith(">"): cur = line[1:].split()[0]; L[cur] = 0
        else: L[cur] += len(line.strip())
    for pid, n in L.items():
        if "-g" not in pid: continue
        fam = pid.split(".", 1)[1].rsplit("-g", 1)[0]
        has = rules[next(k for k in rules if k.lower() == fam)][2]   # the stub hits all its models
        for f, (mlen, rule, need) in rules.items():   # a full ABC transporter is also an ABC_NBD
            if n >= mlen and (need <= has if rule == "ALL" else need & has):
                exp[f] = exp.get(f, 0) + 1
    ok &= all(int(counts[sp][f]) == exp.get(f, 0) for f in rules)
chk("every family count equals the answer computed from the headers and the map", ok)
chk("the dashboard was written", os.path.exists(f"{T}/out/dashboard/index.html"))
chk("one species data file per selected species", len(glob.glob(f"{T}/out/dashboard/data/species/*.js")) == 3)
chk("the release's Pfam version travels with the database", "99.0" in open(f"{T}/out/db/PFAM_VERSION").read())
sys.exit(1 if bad else 0)
PY
grep -q "all steps OK" "$TMP/run.log" && echo "  OK   every step ran"

echo "== --taxdump builds the taxonomy; --no-collapse uses the files as they are =="
mkdir -p "$TMP/taxdump"
printf '1\t|\t1\t|\tno rank\t|\n2\t|\t1\t|\tphylum\t|\n3\t|\t2\t|\tclass\t|\n4\t|\t3\t|\torder\t|\n5\t|\t4\t|\tspecies\t|\n6\t|\t4\t|\tspecies\t|\n' > "$TMP/taxdump/nodes.dmp"
printf '1\t|\troot\t|\t\t|\tscientific name\t|\n2\t|\tArthropoda\t|\t\t|\tscientific name\t|\n3\t|\tInsecta\t|\t\t|\tscientific name\t|\n4\t|\tLepidoptera\t|\t\t|\tscientific name\t|\n5\t|\tGenusins070 species070\t|\t\t|\tscientific name\t|\n5\t|\tseventy moth\t|\t\t|\tgenbank common name\t|\n6\t|\tGenusins071 species071\t|\t\t|\tscientific name\t|\n' > "$TMP/taxdump/names.dmp"
$D run --proteomes "$TMP/fx/proteomes" --out "$TMP/tx" --species Genusins070,Genusins071 --taxdump "$TMP/taxdump" \
   --no-collapse --no-figures --db-url "$TMP/rel/defensome_db.tar.gz" >"$TMP/tx.log" 2>&1 || { tail -30 "$TMP/tx.log"; exit 1; }
python3 - "$TMP/tx" <<'PY'
import sys, re, json, os
T = sys.argv[1]
P = json.loads(re.search(r'const D = (\{.*?\});\n', open(f"{T}/dashboard/index.html").read(), re.S).group(1))
sp = {s["name"]: s for s in P["species"]}
ok = (P["trait_name"] == "class" and sp["Genusins070_species070_pep"]["tax"].get("order") == "Lepidoptera"
      and sp["Genusins070_species070_pep"]["tax"].get("common_name") == "seventy moth"
      and not os.path.exists(f"{T}/samplesheet.tsv"))
print(("  OK   " if ok else "  FAIL ") + "taxonomy and common name from the taxdump; grouped by class; no collapse")
sys.exit(0 if ok else 1)
PY

echo "== gzipped proteomes are read as they come =="
mkdir -p "$TMP/gz"; for f in Genusins070_species070_pep Genusins071_species071_pep; do
  gzip -c "$TMP/fx/proteomes/$f.faa" > "$TMP/gz/$f.fa.gz"; done
$D run --proteomes "$TMP/gz" --out "$TMP/gzo" --no-figures --db-url "$TMP/rel/defensome_db.tar.gz" \
   >"$TMP/gz.log" 2>&1 || { tail -30 "$TMP/gz.log"; exit 1; }
[[ $(tail -n +2 "$TMP/gzo/counts.tsv" | wc -l) -eq 2 ]] && cmp -s <(tail -n +2 "$TMP/tx/counts.tsv") <(tail -n +2 "$TMP/gzo/counts.tsv") \
  && echo "  OK   .fa.gz proteomes give the same counts as plain ones"

echo "== mistakes are caught before any work =="
if $D run --proteomes "$TMP/fx/proteomes" --out "$TMP/o2" --species "Homo sapiens" >"$TMP/r2.log" 2>&1; then
  echo "  FAIL an unknown species was accepted"; exit 1; fi
grep -q "no proteome matches Homo sapiens" "$TMP/r2.log" && echo "  OK   an unknown species is named, and nothing runs"
$D run --proteomes "$TMP/fx/proteomes" --out "$TMP/o3" --species Genusins070 --review \
   --db-url "$TMP/rel/defensome_db.tar.gz" >"$TMP/r3.log" 2>&1
[[ -s "$TMP/o3/samplesheet.tsv" && ! -d "$TMP/o3/hmmsearch" ]] && echo "  OK   --review stops at the sample sheet"
if $D run --proteomes "$TMP/fx/proteomes" --out "$TMP/o4" --species Genusins070 \
     --db-url "$TMP/nowhere.tar.gz" >"$TMP/r4.log" 2>&1; then echo "  FAIL"; exit 1; fi
grep -q "could not download the defensome HMMs" "$TMP/r4.log" && echo "  OK   no database: the error says how to get one"
echo; echo "all run tests passed"
