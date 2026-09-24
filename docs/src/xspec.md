# Using in XSPEC

This page assumes the package is already built ([Building](build.md)) and that
you start `xspec` from the repository root (or `./run-xspec.sh` on older
Linux). Ready-to-run versions of every session below are in
[`examples/xspec/`](https://github.com/phajy/GradusXSPEC.jl/tree/main/examples/xspec).

## Loading the package

```text
XSPEC12> lmod gradusxspec .
```

This loads `libgradusxspec` and registers all models from `model.dat`. The
first Gradus evaluation in a session takes longer while Julia warms up; kernels
computed during a session are cached (see [Performance and caching](@ref)).

## First look at a model

Evaluate a model on a dummy response and plot it
(`examples/xspec/lamp_thin_demo.xcm`):

```text
XSPEC12> dummyrsp 0.5 10.0 300 log
XSPEC12> model gradus_lamp_thin
        (accept defaults with /*)
XSPEC12> newpar 3 2.0        # h = 2 r_g: broader red wing
XSPEC12> setplot energy
XSPEC12> cpd /xs
XSPEC12> plot model
```

Parameter numbering for `gradus_lamp_thin`: `1` spin, `2` inc, `3` h,
`4`–`8` xillver (`Refl_Gamma`, `Refl_A_Fe`, `Refl_logXi`, `Refl_Dens`,
`Refl_Incl`), `9` redshift, `10` norm. The other models follow the same
pattern with their own geometry parameters first (see [Models](models.md)).

For a fast load check that does no ray tracing at all, use the diagnostic
model (`examples/xspec/smoke_test.xcm`):

```text
XSPEC12> dummyrsp 0.3 10.0 200 lin
XSPEC12> model test_gauss
XSPEC12> flux 0.5 9.0
```

## Combining with a continuum

The models are additive components and combine as usual:

```text
XSPEC12> model powerlaw + gradus_lamp_thin
```

A common convention is to tie the reflection photon index to the continuum:

```text
XSPEC12> newpar 6 = 1        # Refl_Gamma tracks powerlaw PhoIndex
```

## Simulating and fitting

`examples/xspec/fake_and_fit.xcm` is a complete template: it freezes the
reflection parameters, simulates a spectrum with `fakeit` at h = 5 r_g,
displaces the height, and refits. The important part for real fits:

```text
XSPEC12> freeze 3-11         # freeze the reflection component ...
XSPEC12> thaw 5              # ... except lamppost height
XSPEC12> fit 100
```

Fitting guidance:

- **Geometry parameters are the expensive ones.** Changing `spin`, `inc`, `h`,
  or `r` moves between corners of the ray-tracing grid; each *new* corner costs
  one kernel computation (about a minute for lamppost, more for ring/disc
  coronae). Revisited corners are served from cache, so fits speed up
  dramatically as the cache fills — and later fits in the same region are
  nearly free once the disk cache is warm.
- **Reflection parameters are cheap.** `Refl_*` only re-interpolate the xillver
  table; vary them freely.
- **Freeze what you do not need.** In particular `inc` and `spin` if they are
  not the target of the analysis.
- Kernels are interpolated multilinearly between grid corners, so best-fit
  values need not lie on the grid.

## Model init string

The trailing field of each model line in `model.dat` is passed to the library
as an init string. The first token is the reflection table path (default
`xillverD-5.fits`, resolved relative to the working directory); additional
tokens enable diagnostics:

| Token | Effect |
|-------|--------|
| `verbose` | Log every evaluation (parameters, cache hits) to the terminal |
| `monitor` | Write fit-progress diagnostics to `gradusxspec_monitor.txt` |
| `monitor=/path/file` | Monitor to a specific file |
| `monitor_interval=N` | Refresh the monitor file every N evaluations (default 10) |

For example, edit the `gradus_lamp_thin` line in `model.dat` to end with
`xillverD-5.fits verbose monitor` and re-run `lmod`.

## Performance and caching

Three cache layers sit between XSPEC and the ray tracer:

1. **L(g) kernels** — RAM plus on-disk persistence (survives XSPEC restarts).
2. **Ring emissivities** — RAM; shared between ring/disc corona evaluations.
3. **Convolution matrices / line spectra** — RAM.

RAM caches share a budget of `GRADUSXSPEC_CACHE_LIMIT_GB` (default 16 GiB).
The disk kernel cache lives in `~/.julia/gradusxspec/kernels` by default and
can be deleted at any time; it is invalidated automatically when the kernel
computation or backend configuration changes.

### Environment variables

Set these before starting `xspec`.

| Variable | Default | Purpose |
|----------|---------|---------|
| `GRADUSXSPEC_VERBOSE` | `0` | Per-evaluation logging (same as the `verbose` init token) |
| `GRADUSXSPEC_CONVOLVE` | `fft` | Blur implementation: `fft` (log-energy FFT) or `matrix` (legacy) |
| `GRADUSXSPEC_FFT_NBINS` | auto | Override the FFT log-grid size |
| `GRADUSXSPEC_BLUR_NATIVE` | `0` | Use the native table grid for blurring instead of the hybrid grid |
| `GRADUSXSPEC_BLUR_EMIN` / `GRADUSXSPEC_BLUR_EMAX` | `2.0` / `150.0` | Core band (keV) of the blur working grid (see [Models — blur grid](models.md#blur-working-grid-and-soft-x-ray-band)) |
| `GRADUSXSPEC_BLUR_DE_ABS` / `GRADUSXSPEC_BLUR_DE_REL` | `0.1` / `0.01` | Bin width caps (keV / fractional) in the core band |
| `GRADUSXSPEC_CACHE_LIMIT_GB` | `16` | Shared RAM budget for all in-memory caches (`0` = unlimited) |
| `GRADUSXSPEC_KERNEL_CACHE` | `1` | Persist L(g) kernels to disk |
| `GRADUSXSPEC_KERNEL_CACHE_DIR` | `~/.julia/gradusxspec/kernels` | Disk kernel cache location |
| `GRADUSXSPEC_MONITOR` | off | Fit monitor (same as the `monitor` init token) |
| `GRADUSXSPEC_MONITOR_INTERVAL` | `10` | Monitor refresh interval (evaluations) |
| `GRADUSXSPEC_KERRZ` | `~/GitHub/kerrz/zig-out/bin/kerrz` | kerrz executable for the `kerrz_*` models |
| `GRADUSXSPEC_KERRZ_NTHREADS` | Julia threads | Threads passed to the kerrz CLI |
| `GRADUSXSPEC_KERRZ_NPHOTONS_LAMP` | `3000` | kerrz lamppost emissivity photon count |
| `GRADUSXSPEC_KERRZ_NPHOTONS_RING` | `50000` | kerrz ring emissivity photon count |
| `GRADUSXSPEC_KERRZ_ROUT` | `400` | kerrz line-profile outer disc radius (r_g) |
| `GRADUSXSPEC_KERRZ_EM_CACHE_DIR` | `~/.julia/gradusxspec/kerrz_em` | kerrz emissivity FITS cache |

## Troubleshooting

- `lmod` fails with a `GLIBCXX` error (Rocky/RHEL 8): start XSPEC with
  `./run-xspec.sh` — see [Building](build.md).
- Segfault in `libtbbmalloc` on macOS: rebuild with the current
  `./build-julia.sh`, which stubs the TBB malloc proxy — see
  [Building](build.md).
- A model evaluates to zero everywhere: check that the reflection table path
  in `model.dat` resolves from the directory where you started `xspec`.
- Blurred flux is zero below ~1.3 keV: the default blur core starts at 2 keV;
  see [Models — blur grid](models.md#blur-working-grid-and-soft-x-ray-band).
  Lower `GRADUSXSPEC_BLUR_EMIN` or set `GRADUSXSPEC_BLUR_NATIVE=1` if you need
  soft-band output.
- To watch what a slow fit is doing, add `monitor` to the init string or set
  `GRADUSXSPEC_MONITOR=1` and follow `gradusxspec_monitor.txt`.
