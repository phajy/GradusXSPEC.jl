# GradusXSPEC.jl

Call Gradus models from XSPEC.

## Models

| XSPEC name | Corona | Disc | Gradus params |
|------------|--------|------|---------------|
| `gradus_lamp_ss` | Lamppost | Shakura–Sunyaev (thick) | spin, Eddington, inc, h + reflection |
| `gradus_lamp_thin` | Lamppost | Thin disc | spin, inc, h + reflection |
| `gradus_ring_thin` | Ring | Thin disc | spin, inc, r, h + reflection |

These models convolve an xillver reflection table with a Gradus relativistic line profile. Load the package once with `lmod gradusxspec .`, then use any model name in `model`.

## Prerequisites

- Julia 1.12 (tested with 1.12.6)
- HEASOFT / XSPEC with `HEADAS` environment configured
- Gradus.jl (installed via this project's `Project.toml`)
- Homebrew `gcc@14` (for XSPEC local-model linking on macOS)

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

### 2. Create the XSPEC model package

```sh
./build-xspec.sh
```

This runs `clean-xspec-package.sh`, `initpackage`, `patch-xspec-makefile.sh`, `fix-heasoft-f77libs.sh`, and `hmake`. Use `./build-xspec.sh --full` for a full clean rebuild (removes the existing `Makefile` and library first). Requires `HEADAS` to be set.

Manual steps, if preferred:

```sh
./clean-xspec-package.sh
initpackage gradusxspec model.dat .
./patch-xspec-makefile.sh
./fix-heasoft-f77libs.sh
hmake
```

`initpackage` generates a fresh `Makefile` for the local model package. If you re-run it after changing `model.dat`, remove stale `lpack_<package>.*` files first with `./clean-xspec-package.sh`. That Makefile does not know about the GradusXSPEC shared library, so `patch-xspec-makefile.sh` inserts the required rpath and link flags into `HD_SHLIB_LIBS`. The script is idempotent (safe to run twice) and anchors on the `-lXS` line rather than fixed line numbers, so it should survive minor HEASOFT/Makefile changes better than a static patch file.

After upgrading Homebrew `gcc@14`, `./build-xspec.sh` also refreshes HEASOFT's Fortran library paths automatically via `fix-heasoft-f77libs.sh`.

### 3. Test the models

```
xspec
lmod gradusxspec .
model gradus_lamp_ss
```

Set the model parameters, evaluate over a dummy energy grid, and use `iplot` to visualize the results. Use `model gradus_lamp_thin` for the thin-disc lamppost variant (no `Eddington` parameter), or `model gradus_ring_thin` for a co-rotating ring corona (`spin`, `inc`, `r`, `h`).

```
========================================================================
Model gradus_lamp_ss<1> Source No.: 1   Active/Off
Model Model Component  Parameter  Unit     Value
 par  comp
   1    1   gradus_lamp_ss spin                0.900000     +/-  0.0
   2    1   gradus_lamp_ss Eddington           0.150000     +/-  0.0
   3    1   gradus_lamp_ss inc        degrees  40.0000      +/-  0.0
   4    1   gradus_lamp_ss h          r_g      6.00000      +/-  0.0
   ...
```

![XSPEC plot output](figs/pgplot.png)

## Diagnostics

Terminal logging and a fit-monitor file can be enabled via environment variables or tokens in the XSPEC model init string (after the table path).

### Environment variables

| Variable | Effect |
|----------|--------|
| `GRADUSXSPEC_VERBOSE=1` | Print per-evaluation parameter and cache messages to the terminal (`true` / `yes` also work) |
| `GRADUSXSPEC_MONITOR=1` | Write fit diagnostics to `gradusxspec_monitor.txt` in the repo root |
| `GRADUSXSPEC_MONITOR=/path/to/file` | Same, but use a custom output path |
| `GRADUSXSPEC_MONITOR_INTERVAL=N` | Refresh the monitor file every `N` evaluations (default 10) |

Example:

```sh
export GRADUSXSPEC_VERBOSE=1
export GRADUSXSPEC_MONITOR=1
xspec
```

### Init-string tokens

The model init string starts with the reflection table path (default `xillverD-5.fits`). Extra tokens may follow:

| Token | Effect |
|-------|--------|
| `verbose` | Same as `GRADUSXSPEC_VERBOSE` |
| `monitor` | Enable the monitor file at the default path |
| `monitor=/path/to/file` | Enable the monitor file at a custom path |
| `monitor_interval=N` | Refresh interval (default 10) |

In XSPEC, set these on the model’s init string, for example:

```
xillverD-5.fits verbose monitor
```

The monitor file records current parameters, interpolation brackets, fit-history histograms, and cache hit rates. It is gitignored.

## References

- [Appendix C: Adding Models to XSPEC](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSappendixLocal.html) — official HEASARC documentation for local models (`initpackage`, `model.dat`, `hmake`, `lmod`)
- [OGIP 92-009: The XSPEC Table Model Format (PDF)](https://heasarc.gsfc.nasa.gov/FTP/caldb/docs/memos/ogip_92_009/ogip_92_009.pdf) — reference format for table model files, relevant when using table-file contents in model evaluation
- [Gradus.jl codebase](https://codeberg.org/astro-group/Gradus.jl) — upstream source used to build model physics and APIs
- [Gradus.jl documentation](https://astro-group.codeberg.page/Gradus.jl/dev/) — reference for Gradus usage, line profiles, emissivity, and API behavior
