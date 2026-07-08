#!/usr/bin/env bash
# Insert GradusXSPEC rpath/link flags into the XSPEC initpackage Makefile.
# Idempotent: safe to run more than once.
set -euo pipefail

cd "$(dirname "$0")"

MAKEFILE="${1:-Makefile}"
MARKER="-lGradusXSPEC"

if [[ ! -f "$MAKEFILE" ]]; then
  echo "error: ${MAKEFILE} not found (run 'initpackage gradus model.dat .' first)" >&2
  exit 1
fi

if grep -qF -- "$MARKER" "$MAKEFILE"; then
  echo "${MAKEFILE} already configured for GradusXSPEC"
  exit 0
fi

python3 - "$MAKEFILE" <<'PY'
import re
import sys

makefile = sys.argv[1]
marker = "-lGradusXSPEC"
flags = [
    "-Wl,-rpath,${PWD}/build/lib",
    "-Wl,-rpath,${PWD}/build/lib/julia",
    "-L${PWD}/build/lib",
    marker,
]

with open(makefile) as f:
    lines = f.readlines()

if any(marker in line for line in lines):
    print(f"{makefile} already configured for GradusXSPEC")
    sys.exit(0)

out = []
patched = False

for i, line in enumerate(lines):
    out.append(line)

    if patched:
        continue

    if "-lXS" in line and line.rstrip().endswith("\\"):
        indent = "\t\t\t\t\t\t  "
        if i + 1 < len(lines):
            m = re.match(r"^[ \t]*", lines[i + 1])
            if m and m.group(0):
                indent = m.group(0)
        for flag in flags:
            out.append(f"{indent}{flag} \\\n")
        patched = True

if not patched:
    print(
        "error: could not find HD_SHLIB_LIBS continuation line containing '-lXS'",
        file=sys.stderr,
    )
    sys.exit(1)

with open(makefile, "w") as f:
    f.writelines(out)

print(f"Configured {makefile} for GradusXSPEC")
PY
