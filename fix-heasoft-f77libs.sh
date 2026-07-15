#!/usr/bin/env bash
# Refresh HEASOFT hmakerc F77LIBS4C paths for the installed Homebrew gcc@14.
#
# macOS-only workaround. HEASOFT records absolute Fortran library paths at
# configure time. After `brew upgrade gcc@14`, local model links can fail with
# "library emutls_w not found". This script rewrites F77LIBS4C in every hmakerc
# under the HEASOFT tree using paths reported by gfortran-14.
#
# On other platforms (or without Homebrew) this is a no-op.
#
# Usage (HEADAS must be set):
#   ./fix-heasoft-f77libs.sh
#
# Optional:
#   HEASOFT_ROOT=/path/to/heasoft-6.36 ./fix-heasoft-f77libs.sh

set -euo pipefail

cd "$(dirname "$0")"

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "fix-heasoft-f77libs: not macOS; nothing to do"
  exit 0
fi

if [[ -z "${HEADAS:-}" ]]; then
  echo "error: HEADAS is not set (source \$HEADAS/headas-init.sh)" >&2
  exit 1
fi

if ! command -v brew >/dev/null 2>&1; then
  echo "fix-heasoft-f77libs: Homebrew not found; nothing to do"
  exit 0
fi

GCC14_PREFIX="$(brew --prefix gcc@14 2>/dev/null || true)"
if [[ -z "$GCC14_PREFIX" || ! -d "$GCC14_PREFIX" ]]; then
  echo "error: gcc@14 is not installed (brew install gcc@14)" >&2
  exit 1
fi

GFORTRAN="${GCC14_PREFIX}/bin/gfortran-14"
if [[ ! -x "$GFORTRAN" ]]; then
  GFORTRAN="$(command -v gfortran-14 || true)"
fi
if [[ -z "$GFORTRAN" || ! -x "$GFORTRAN" ]]; then
  echo "error: gfortran-14 not found (expected under ${GCC14_PREFIX}/bin)" >&2
  exit 1
fi

EMUTLS_LIB="$("$GFORTRAN" -print-file-name=libemutls_w.a)"
GFORTRAN_LIB="$("$GFORTRAN" -print-file-name=libgfortran.dylib)"
GCC_LIB="$("$GFORTRAN" -print-file-name=libgcc.a)"

for lib in "$EMUTLS_LIB" "$GFORTRAN_LIB" "$GCC_LIB"; do
  if [[ ! -f "$lib" ]]; then
    echo "error: expected library not found: $lib" >&2
    exit 1
  fi
done

EMUTLS_DIR="$(cd "$(dirname "$EMUTLS_LIB")" && pwd)"
GCC_LIB_DIR="$(cd "$(dirname "$EMUTLS_DIR")" && pwd)"
GFORTRAN_DIR="$(cd "$(dirname "$GFORTRAN_LIB")" && pwd)"

F77LIBS4C="-L/usr/lib -L${EMUTLS_DIR} -L${GCC_LIB_DIR} -L${GFORTRAN_DIR} -lemutls_w -lheapt_w -lgfortran -lquadmath "

HEASOFT_ROOT="${HEASOFT_ROOT:-$(dirname "$HEADAS")}"
if [[ ! -d "$HEASOFT_ROOT" ]]; then
  echo "error: HEASOFT root not found: $HEASOFT_ROOT" >&2
  exit 1
fi

echo "HEASOFT root: $HEASOFT_ROOT"
echo "gfortran:     $GFORTRAN"
echo "F77LIBS4C:    $F77LIBS4C"

updated=0
skipped=0
found=0
expected_line="F77LIBS4C=\"${F77LIBS4C}\""

while IFS= read -r file; do
  [[ -z "$file" ]] && continue
  found=$((found + 1))
  current_line="$(grep '^F77LIBS4C=' "$file" || true)"
  if [[ "$current_line" == "$expected_line" ]]; then
    skipped=$((skipped + 1))
    continue
  fi
  python3 - "$file" "$F77LIBS4C" <<'PY'
import re
import sys

path, value = sys.argv[1], sys.argv[2]
with open(path) as f:
    text = f.read()

new_text, n = re.subn(
    r'^F77LIBS4C=".*"$',
    f'F77LIBS4C="{value}"',
    text,
    count=1,
    flags=re.MULTILINE,
)
if n != 1:
    print(f"error: could not update F77LIBS4C in {path}", file=sys.stderr)
    sys.exit(1)

with open(path, "w") as f:
    f.write(new_text)
PY
  updated=$((updated + 1))
done < <(find "$HEASOFT_ROOT" -name hmakerc -exec grep -l '^F77LIBS4C=' {} + 2>/dev/null || true)

if [[ $found -eq 0 ]]; then
  echo "error: no hmakerc files with F77LIBS4C found under $HEASOFT_ROOT" >&2
  exit 1
fi

echo "hmakerc files: ${found} total, ${updated} updated, ${skipped} already current"
