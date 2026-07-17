#!/usr/bin/env bash
# Build the XSPEC local model package for GradusXSPEC.
#
# Chains: clean stale initpackage files -> initpackage -> patch Makefile ->
#         fix HEASOFT F77 libs -> hmake
#
# Prerequisites:
#   - HEADAS set (e.g. source $HEADAS/headas-init.sh)
#   - ./build-julia.sh already run (produces model.dat and libGradusXSPEC)
#
# Usage (from repository root):
#   ./build-xspec.sh
#   ./build-xspec.sh --full   # full clean before initpackage (Makefile, library, etc.)

set -euo pipefail

cd "$(dirname "$0")"

FULL_CLEAN=false
if [[ "${1:-}" == "--full" ]]; then
  FULL_CLEAN=true
elif [[ -n "${1:-}" ]]; then
  echo "usage: $0 [--full]" >&2
  exit 1
fi

if [[ -z "${HEADAS:-}" ]]; then
  echo "error: HEADAS is not set (source \$HEADAS/headas-init.sh)" >&2
  exit 1
fi

if [[ ! -f model.dat ]]; then
  echo "error: model.dat not found (run ./build-julia.sh first)" >&2
  exit 1
fi

PACKAGE_NAME="$(julia --project=. -e 'include("src/model_definition.jl"); print(PACKAGE_NAME)')"
MODEL_DAT="model.dat"

for cmd in initpackage hmake; do
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "error: '$cmd' not found on PATH (is HEADAS initialized?)" >&2
    exit 1
  fi
done

echo "==> Cleaning stale initpackage artifacts"
if [[ "$FULL_CLEAN" == true ]]; then
  ./clean-xspec-package.sh --full
else
  ./clean-xspec-package.sh
fi

echo "==> Running initpackage ${PACKAGE_NAME} ${MODEL_DAT} ."
initpackage "${PACKAGE_NAME}" "${MODEL_DAT}" .

echo "==> Patching Makefile for GradusXSPEC"
./patch-xspec-makefile.sh

# macOS-only workaround: Homebrew gcc@14 upgrades break HEASOFT's recorded
# Fortran library paths. Rather than mutating the shared HEASOFT install, we
# compute the correct paths and pass them to hmake on the command line, which
# takes precedence over the value baked into hmakerc and scopes the fix to this
# build. Linux HEASOFT installs do not need (or have) brew, so this is skipped.
HMAKE_ARGS=()
if [[ "$(uname -s)" == "Darwin" ]] && command -v brew >/dev/null 2>&1; then
  echo "==> Computing HEASOFT Fortran library paths (gcc@14)"
  F77LIBS4C_OVERRIDE="$(./fix-heasoft-f77libs.sh --print)"
  if [[ -n "$F77LIBS4C_OVERRIDE" ]]; then
    HMAKE_ARGS+=("F77LIBS4C=${F77LIBS4C_OVERRIDE}")
    echo "    F77LIBS4C=${F77LIBS4C_OVERRIDE}"
  fi
else
  echo "==> Skipping HEASOFT gcc@14 path fix (not macOS with Homebrew)"
fi

echo "==> Building with hmake"
# Guard the array expansion for bash 3.2 (stock macOS) under `set -u`.
hmake ${HMAKE_ARGS[@]+"${HMAKE_ARGS[@]}"}

echo "Done. Load in XSPEC with: lmod ${PACKAGE_NAME} ."
echo "Models: gradus_lamp_ss, gradus_lamp_thin, gradus_ring_thin, gradus_disc_thin, test_gauss"
