#!/usr/bin/env bash
# Headless test of the dashboard JavaScript against a real generated dashboard.
# Needs node. Usage:  bash test/test_dashboard.sh results/dashboard/index.html
set -euo pipefail
HTML=${1:?usage: bash test/test_dashboard.sh <dashboard/index.html>}
HERE=$(cd "$(dirname "$0")" && pwd)
command -v node >/dev/null || { echo "node not found; skipping"; exit 0; }
TMP=$(mktemp -d)
python3 - "$HTML" "$TMP" <<'PY'
import re, sys, json
html = open(sys.argv[1]).read()
m = re.search(r'const D = (\{.*?\});\n', html, re.S)
if not m: sys.exit("could not find the embedded payload")
json.loads(m.group(1))
open(sys.argv[2] + "/payload.js", "w").write("global.D = " + m.group(1) + ";\n")
i = html.index("const D = "); j = html.index("</script>", i)
k = html.index("<script>", j); e = html.index("</script>", k)
open(sys.argv[2] + "/app.js", "w").write(html[k+8:e])
print("payload parses; JS extracted")
PY
cat "$HERE/dom_shim.js" "$TMP/payload.js" "$TMP/app.js" "$HERE/dashboard_tests.js" \
    "$HERE/tree_tests.js" > "$TMP/run.js"
node "$TMP/run.js"
rm -rf "$TMP"
