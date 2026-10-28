#!/usr/bin/env bash
# Taxonomic presentation: name display, label decluttering, habitat sort/filter.
# Needs node + playwright; without them the fixture is still built and the
# payload checked, which catches the metadata column selection on its own.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
python3 "$HERE/make_taxonomy_fixture.py" "$TMP/res"
LOG=$TMP/dash.log
python3 "$ROOT/defensome.py" dashboard --out "$TMP/res" \
  --metadata "$TMP/res/metadata.tsv" --group-by class >"$LOG" 2>&1 || { cat "$LOG"; exit 1; }

# habitat must survive column selection; common_name must not become a filter
grep -q "taxonomy columns: phylum, class, order, habitat" "$LOG" || {
  echo "  FAIL metadata columns:"; grep -i "taxonomy columns" "$LOG"; exit 1; }
echo "  OK   habitat is picked up as a grouping column, common_name is not"

bash "$HERE/test_dashboard.sh" "$TMP/res/dashboard/index.html"
if node -e "require(process.env.PLAYWRIGHT_MODULE || 'playwright')" 2>/dev/null; then
  echo; echo "== real browser =="
  node "$HERE/taxonomy_browser_test.js" "$TMP/res/dashboard/index.html"
  echo; echo "== layout at catalogue scale (1,236 species, 32 classes) =="
  python3 "$HERE/make_taxonomy_fixture.py" "$TMP/cat" --catalogue >/dev/null
  python3 "$ROOT/defensome.py" dashboard --out "$TMP/cat" \
    --metadata "$TMP/cat/metadata.tsv" --group-by class >/dev/null 2>&1
  node "$HERE/layout_browser_test.js" "$TMP/cat/dashboard/index.html"
else
  echo; echo "playwright not installed; skipping the real-browser tests"
  echo "  npm i playwright && npx playwright install chromium"
fi
