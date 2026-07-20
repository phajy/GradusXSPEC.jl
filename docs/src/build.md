# Building

## Prerequisites

- Julia 1.12 (tested with 1.12.6)
- HEASOFT / XSPEC with `HEADAS` configured
- Gradus.jl (via this project's `Project.toml`)
- `xillverD-5.fits` in the repository root ([download](https://sites.srl.caltech.edu/~javier/xillver/index.html))
- Homebrew `gcc@14` on **macOS only** (Fortran linking workaround)

## Check the environment

Before building, run the read-only diagnostic:

```sh
./check-env.sh
```

It reports PASS / WARN / FAIL for Julia, `HEADAS`, HEASOFT tools, build inputs, and
(on macOS) the gcc@14 Fortran libraries.

## Build steps

### 1. Julia shared library

From the repository root:

```sh
./build-julia.sh
```

This instantiates the Julia environment, compiles `libGradusXSPEC`, and writes
`model.dat`. The library is at `build/lib/libGradusXSPEC.{dylib,so}`.

### 2. XSPEC local model package

With `HEADAS` sourced:

```sh
./build-xspec.sh
```

This runs `clean-xspec-package.sh`, `initpackage`, `patch-xspec-makefile.sh`, and
`hmake`. Use `./build-xspec.sh --full` for a full clean rebuild.

`patch-xspec-makefile.sh` inserts rpath and link flags for the GradusXSPEC
library into the generated Makefile (idempotent, anchored on the `-lXS` line).

### 3. Load in XSPEC

```sh
./run-xspec.sh   # preferred on older Linux; plain `xspec` is fine elsewhere
```

Then in the XSPEC prompt:

```
lmod gradusxspec .
model gradus_lamp_ss
```

Evaluate over a dummy response and inspect with `iplot`. On Rocky/RHEL 8, use
`./run-xspec.sh` if `lmod` fails with a `GLIBCXX` error (see below).

## Platform notes

### Linux (native HEASoft)

On Linux the build pipeline is the same as macOS except the gcc@14 workaround is
skipped entirely — a standard HEASOFT install's Fortran paths are used as-is.

```sh
source $HEADAS/headas-init.sh
./check-env.sh
./build-julia.sh
./build-xspec.sh
```

This is the recommended path if you already have HEASOFT installed.

#### Julia `libstdc++` on older Linux (Rocky / RHEL 8)

Julia 1.12 needs a newer `libstdc++` than the system library on Rocky/RHEL 8
(`/lib64/libstdc++.so.6` typically tops out at `GLIBCXX_3.4.25`). Loading
`libgradusxspec.so` can then fail with:

```text
/lib64/libstdc++.so.6: version `GLIBCXX_...' not found
(required by .../build/lib/julia/libjulia-internal.so...)
```

juliaup’s Julia install ships a suitable `libstdc++` under its `lib/` directory,
but XSPEC does not inherit that path: HEASoft’s `lib/` and `/lib64` are searched
first. Use the wrapper so Julia’s libraries are prepended to `LD_LIBRARY_PATH`
before starting XSPEC:

```sh
./run-xspec.sh
# or a smoke script:
./run-xspec.sh - docker/smoke-test.xcm
```

Then `lmod gradusxspec .` as usual. On newer distributions whose system
`libstdc++` already provides the needed GLIBCXX versions, plain `xspec` is fine.

### macOS (gcc@14 Fortran linking)

HEASOFT records absolute Fortran library paths at configure time. After
`brew upgrade gcc@14`, local-model links can fail with `library emutls_w not
found`. `build-xspec.sh` computes the correct paths with
`fix-heasoft-f77libs.sh --print` and passes them to `hmake` on the command line,
scoping the fix to this build without mutating the shared HEASOFT install.

As a last resort, `./fix-heasoft-f77libs.sh` (no arguments) rewrites `F77LIBS4C`
across the HEASOFT tree in place.

### Docker (optional reproducibility)

If you do not have a local Linux HEASOFT install, or want a clean-room check, see
[`docker/README.md`](https://github.com/phajy/GradusXSPEC.jl/blob/main/docker/README.md).
Docker is optional; native Linux builds are the usual development path.

## Building the manual

```sh
./build-docs.sh
```

Open `docs/build/index.html` in a browser.

## Continuous integration

A practical split for GitHub Actions:

| Tier | What | When | Cost |
|------|------|------|------|
| **Light** | `Pkg.instantiate()`, load `GradusXSPEC`, Julia validation scripts | Every push / PR | Minutes |
| **Heavy** | Full HEASoft + PackageCompiler + `hmake` + XSPEC smoke test | Manual (`workflow_dispatch`) or release tags | 30–60+ min, ~12 GB image |

A full end-to-end XSPEC test in CI is **feasible but expensive**: the HEASoft
Docker image is large, PackageCompiler is slow, and runners need enough disk and
memory. Running it rarely (weekly, pre-release, or on demand) is reasonable;
running it on every commit is probably not worth the queue time and storage on
GitHub-hosted runners unless you use a self-hosted Linux machine with HEASOFT
already installed.

Recommended starting point:

1. **Julia-only CI** on every PR — `validate_table.jl`, `validate_spectrum.jl`
   (with a cached or downloaded `xillverD-5.fits`), no XSPEC.
2. **Optional HEASoft workflow** triggered manually — build the package and run
   `docker/smoke-test.xcm` or the equivalent native steps on a self-hosted runner.

## Diagnostics

Environment variables and init-string tokens for verbose logging and fit
monitoring are documented in the repository [README](https://github.com/phajy/GradusXSPEC.jl#diagnostics).
