#!/usr/bin/env bash
# Build the release tarball with correct file modes.
# cp and tar do not reliably carry the executable bit, and the test suite and
# the SLURM scripts both depend on it.
set -euo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$HERE"
chmod +x run_pipeline.sh package.sh slurm/*.sh test/*.sh test/*.py 2>/dev/null || true
rm -rf test/out
python3 - <<'PY'
import base64, zlib, re
s = open("defensome.py").read()
blobs = {n: base64.b64encode(zlib.compress(open(n, "rb").read(), 9)).decode()
         for n in ("defensome_map.tsv", "dashboard.html", "dashboard.js")}
lines = ["_ASSETS = {"]
for n, b in blobs.items():
    lines.append(f'    "{n}": (')
    for i in range(0, len(b), 100):
        lines.append(f'        "{b[i:i+100]}"')
    lines.append("    ),")
lines.append("}")
s = re.sub(r"_ASSETS = \{.*?\n\}\n", "\n".join(lines) + "\n", s, count=1, flags=re.S)
open("defensome.py", "w").write(s)
print("embedded assets refreshed")
PY
python3 -c "import ast; ast.parse(open('defensome.py').read())" && echo "defensome.py parses"
for f in run_pipeline.sh slurm/*.sh test/*.sh; do bash -n "$f"; done && echo "shell scripts parse"
echo "--- test suites ---"
bash test/test_pipeline.sh       > /dev/null && echo "  pipeline OK"
bash test/test_small_proteins.sh > /dev/null && echo "  small proteins OK"
bash test/test_run.sh           | grep -E "passed|FAIL" | sed 's/^/  /'
if node -e "require(process.env.PLAYWRIGHT_MODULE || 'playwright')" 2>/dev/null; then
  bash test/test_trees.sh       | grep -E "passed|FAIL" | sed 's/^/  /'
  bash test/test_taxonomy_ui.sh | grep -E "passed|FAIL" | sed 's/^/  /'
  bash test/test_site.sh        | grep -E "passed|FAIL" | sed 's/^/  /'
else
  echo "  trees and taxonomy UI: node/playwright absent, skipped"
fi
echo "--- version ---"; python3 defensome.py --version
