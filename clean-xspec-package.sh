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
  MODEL_NAME="$(awk 'NF { print $1; exit }' model.dat)"
elif [[ -f src/model_definition.jl ]]; then
  MODEL_NAME="$(julia --project=. -e 'include("src/model_definition.jl"); print(MODEL_NAME)')"
else
  echo "error: could not determine model name (run ./build-julia.sh first)" >&2
  exit 1
fi

echo "Cleaning initpackage artifacts for model: ${MODEL_NAME}"

rm -f \
  "lpack_${MODEL_NAME}.cxx" \
  "lpack_${MODEL_NAME}.o" \
  "${MODEL_NAME}FunctionMap.cxx" \
  "${MODEL_NAME}FunctionMap.h"

if [[ "$FULL_CLEAN" == true ]]; then
  shopt -s nullglob
  rm -f "lib${MODEL_NAME}".* Makefile pkgIndex.tcl *.bck
  shopt -u nullglob
  echo "Full clean: removed Makefile, pkgIndex.tcl, and lib${MODEL_NAME}.*"
fi

echo "Done."
