#!/usr/bin/env bash
# Tree tab: known-answer tests in node, then real-browser interaction tests.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd); ROOT=$(cd "$HERE/.." && pwd)
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT
python3 "$HERE/make_tree_fixture.py" "$TMP/res"
python3 "$ROOT/defensome.py" dashboard --out "$TMP/res" --samplesheet "$TMP/res/samplesheet.tsv" >/dev/null
bash "$HERE/test_dashboard.sh" "$TMP/res/dashboard/index.html"
if node -e "require(process.env.PLAYWRIGHT_MODULE || 'playwright')" 2>/dev/null; then
  echo; echo "== real browser =="; node "$HERE/browser_test.js" "$TMP/res/dashboard/index.html"
else
  echo; echo "playwright not installed; skipping the real-browser tests"
  echo "  npm i playwright && npx playwright install chromium"
fi
