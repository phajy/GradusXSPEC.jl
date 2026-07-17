# Validation

Check the Julia-side pipeline before or alongside XSPEC testing. The scripts
exercise table interpolation, line profiles, convolution, and the full blurred
spectrum.

## Prerequisites

- Repository cloned and `julia --project=.` instantiated
- `xillverD-5.fits` in the repository root
- **Plots.jl** for figures (one-time):

```sh
julia --project=. -e 'using Pkg; Pkg.add("Plots")'
```

## Numeric checks

```sh
julia --project=. scripts/validate_table.jl
julia --project=. scripts/validate_line_profile.jl
julia --project=. scripts/validate_convolution.jl
julia --project=. scripts/validate_spectrum.jl
```

The line-profile script calls Gradus and can take one to two minutes on first run.
`validate_spectrum.jl` compares the direct evaluation path against the
grid-interpolated path used in XSPEC.

## Generate plots

```sh
julia --project=. scripts/plot_validation.jl
```

Figures are written to `docs/figures/` (gitignored). Rebuild this manual with
`./build-docs.sh` to view them inline after generation.

!!! note "Flux units"
    XSPEC table models store **integrated flux per energy bin**. Internal
    calculations use those bin-integrated values. Plots divide by bin width to
    show **flux density (per keV)** for readability on irregular energy grids.

## Default test parameters

| Component | Parameters |
|-----------|------------|
| Gradus (lamppost) | spin=0.998, Eddington=0.1, inc=30°, h=3 r_g |
| xillver reflection | Γ=2.0, A_Fe=1.0, log ξ=2.0, n_e=17, incl=45° |

XSPEC adds `norm` and applies `*redshift` from `model.dat`; the Julia API uses
the nine physics parameters only.

## Expected figures

After running `plot_validation.jl`, expect:

1. **Line-profile kernel `L(g)`** — unit-area kernel peaked near `g ≈ 1`.
2. **Reflection spectrum `R(E)`** — interpolated xillver on the table grid.
3. **Fe Kα region (3–9 keV)** — zoom on the iron-line band.
4. **Narrow-kernel convolution** — near-identity blur; curves should overlap.
5. **Gradus-kernel convolution** — Fe features broadened and red-wing extended.
6. **Flux-conserving rebin** — integrated flux per bin conserved.

## Pipeline

```text
xillver FITS  →  R(E)              [table_model.jl]
Gradus params  →  L(g)              [line_profile.jl]
R, L           →  F = M * R         [convolution.jl]
F              →  XSPEC energies    [spectrum.jl]
```

`evaluate_spectrum_interpolated` is what the XSPEC entry points call; it
interpolates across parameter grids and caches convolution matrices.
`evaluate_spectrum` evaluates at exact parameters and is used for validation.

## XSPEC smoke test (Linux or macOS)

After `./build-xspec.sh`:

```
xspec
lmod gradusxspec .
dummyrsp 0.3 10.0 200 lin
model test_gauss
/*
flux 0.5 9.0
```

The `test_gauss` model is fast (no ray tracing). Physical models use the same
interface but each evaluation takes longer.

!!! note "Documentation status"
    Ring and disc corona validation plots will be added after the Linux build is
    verified on native HEASOFT.
