# GradusXSPEC.jl

Call Gradus models from XSPEC.

## Prerequisites

- Julia 1.12 (tested with 1.12.6)
- HEASOFT / XSPEC with `HEADAS` environment configured
- Gradus.jl (installed via this project's `Project.toml`)

## How to build the library and import it into XSPEC

### 1. Build the Julia library

From the repository root:

```sh
./build-julia.sh
```

This instantiates the Julia environment, then compiles the shared library and writes `model.dat`. The library is placed at `build/lib/libGradusXSPEC.dylib` on macOS (`.so` on Linux).

Alternatively, run the steps manually:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. src/build_lib.jl
```

### 2. Create the XSPEC model

```sh
initpackage gradus model.dat .
./patch-xspec-makefile.sh
hmake
```

`initpackage` generates a fresh `Makefile` for the local model. That Makefile does not know about the GradusXSPEC shared library, so `patch-xspec-makefile.sh` inserts the required rpath and link flags into `HD_SHLIB_LIBS`. The script is idempotent (safe to run twice) and anchors on the `-lXS` line rather than fixed line numbers, so it should survive minor HEASOFT/Makefile changes better than a static patch file.

### 3. Test the model

```
xspec
lmod gradus .
```

Set the model parameters, evaluate over a dummy energy grid, and use `iplot` to visualize the results.

```
========================================================================
Model gradus<1> Source No.: 1   Active/Off
Model Model Component  Parameter  Unit     Value
 par  comp
   1    1   gradus     spin                0.900000     +/-  0.0
   2    1   gradus     Eddington           0.150000     +/-  0.0
   3    1   gradus     inc        degrees  40.0000      +/-  0.0
   4    1   gradus     h          r_g      6.00000      +/-  0.0
   5    1   gradus     norm                1.00000      +/-  0.0
________________________________________________________________________
```

![XSPEC plot output](figs/pgplot.png)
