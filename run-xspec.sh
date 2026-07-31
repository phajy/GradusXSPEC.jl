#!/usr/bin/env bash
# Launch XSPEC with Julia's shared libraries ahead of HEASoft/system paths.
#
# On older Linux distributions (e.g. Rocky/RHEL 8), the system libstdc++ is too
# old for Julia 1.12 (needs GLIBCXX_3.4.26+). Loading libgradusxspec.so then
# fails with:
#   /lib64/libstdc++.so.6: version `GLIBCXX_...' not found
#   (required by .../build/lib/julia/libjulia-internal.so...)
#
# Julia (juliaup) ships a newer libstdc++ under its install tree. This wrapper
# prepends that directory to LD_LIBRARY_PATH after HEADAS is initialized, so
# XSPEC's lmod finds Julia's C++ runtime before /lib64 or HEASoft's lib/.
#
# Usage (from repository root; HEADAS must already be set, or pass headas-init):
#   ./run-xspec.sh
#   ./run-xspec.sh - docker/smoke-test.xcm
#
# Optional:
#   HEADAS_INIT=/path/to/headas-init.sh ./run-xspec.sh

set -euo pipefail

cd "$(dirname "$0")"

if [[ -z "${HEADAS:-}" ]]; then
  if [[ -n "${HEADAS_INIT:-}" && -f "${HEADAS_INIT}" ]]; then
    # shellcheck disable=SC1090
    source "${HEADAS_INIT}"
  else
    echo "error: HEADAS is not set (source \$HEADAS/headas-init.sh, or set HEADAS_INIT)" >&2
    exit 1
  fi
fi

if ! command -v xspec >/dev/null 2>&1; then
  echo "error: xspec not found on PATH (is HEADAS initialized?)" >&2
  exit 1
fi

if ! command -v julia >/dev/null 2>&1; then
  echo "error: julia not found on PATH" >&2
  exit 1
fi

JULIA_ROOT="$(julia -e 'print(dirname(Sys.BINDIR))')"
JULIA_LIB="${JULIA_ROOT}/lib"
JULIA_LIB_JULIA="${JULIA_ROOT}/lib/julia"

for d in "$JULIA_LIB" "$JULIA_LIB_JULIA"; do
  if [[ ! -d "$d" ]]; then
    echo "error: expected Julia library directory missing: $d" >&2
    exit 1
  fi
done

# Prepend after HEADAS so Julia wins over HEASoft's and the system's libstdc++.
export LD_LIBRARY_PATH="${JULIA_LIB}:${JULIA_LIB_JULIA}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

exec xspec "$@"
