#!/usr/bin/env bash
# Build the Documenter manual locally.
#
# Usage (from repository root):
#   ./build-docs.sh
#
# Output: docs/build/index.html

set -euo pipefail

cd "$(dirname "$0")"

julia --project=docs -e 'using Pkg; Pkg.instantiate()'
julia --project=docs docs/make.jl

echo "Done. Open docs/build/index.html"
