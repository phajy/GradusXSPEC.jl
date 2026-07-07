#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. src/build_lib.jl
