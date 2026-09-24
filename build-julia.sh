#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

julia --project=. -e '
using Pkg
Pkg.Registry.add("General")
Pkg.Registry.add(url = "https://github.com/astro-group-bristol/AstroRegistry")
Pkg.instantiate()
'
julia --project=. src/build_lib.jl
