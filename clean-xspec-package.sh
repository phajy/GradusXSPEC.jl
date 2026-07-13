#!/usr/bin/env bash
# Remove XSPEC initpackage artifacts that block re-running initpackage.
#
# initpackage refuses to overwrite files such as lpack_<model>.cxx if they
# already exist. The model name is read from model.dat (first field) when
# present, otherwise from src/model_definition.jl.
#
# Usage (from repository root):
#   ./clean-xspec-package.sh
#   ./clean-xspec-package.sh --full   # also remove Makefile, pkgIndex.tcl, library

set -euo pipefail

cd "$(dirname "$0")"

FULL_CLEAN=false
if [[ "${1:-}" == "--full" ]]; then
  FULL_CLEAN=true
fi

if [[ -f model.dat ]]; then
  PACKAGE_NAME="$(julia --project=. -e 'include("src/model_definition.jl"); print(PACKAGE_NAME)')"
elif [[ -f src/model_definition.jl ]]; then
  PACKAGE_NAME="$(julia --project=. -e 'include("src/model_definition.jl"); print(PACKAGE_NAME)')"
else
  echo "error: could not determine package name (run ./build-julia.sh first)" >&2
  exit 1
fi

echo "Cleaning initpackage artifacts for package: ${PACKAGE_NAME}"

rm -f \
  "lpack_${PACKAGE_NAME}.cxx" \
  "lpack_${PACKAGE_NAME}.o" \
  "${PACKAGE_NAME}FunctionMap.cxx" \
  "${PACKAGE_NAME}FunctionMap.h"

if [[ "$FULL_CLEAN" == true ]]; then
  shopt -s nullglob
  rm -f "lib${PACKAGE_NAME}".* Makefile pkgIndex.tcl *.bck
  # Remove legacy single-model package artifacts from before the rename.
  rm -f lpack_gradus.* gradusFunctionMap.* libgradus.*
  shopt -u nullglob
  echo "Full clean: removed Makefile, pkgIndex.tcl, and lib${PACKAGE_NAME}.*"
fi

echo "Done."
