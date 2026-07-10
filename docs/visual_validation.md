# Visual validation (development)

This page documents how to check the Julia-side model components before wiring
everything into XSPEC. It is aimed at developers working on the reflection +
convolution pipeline.

## Prerequisites

- Repository cloned and `julia --project=.` instantiated
- `xillverD-5.fits` in the repository root ([xillver download](https://sites.srl.caltech.edu/~javier/xillver/index.html))
- **Plots.jl** for figures (one-time):

```sh
julia --project=. -e 'using Pkg; Pkg.add("Plots")'
```

## Quick numeric checks

```sh
julia --project=. scripts/validate_table.jl
julia --project=. scripts/validate_line_profile.jl
julia --project=. scripts/validate_convolution.jl
julia --project=. scripts/validate_spectrum.jl
```

The line-profile script calls Gradus and can take ~1–2 minutes on first run.

## Generate plots

```sh
julia --project=. scripts/plot_validation.jl
```

Figures are written to `docs/figures/`. Open the PNG files or view them inline
below after running the script.

**Note:** XSPEC table models store **integrated flux per energy bin**. All
internal calculations (convolution, rebin, future XSPEC interface) use those
bin-integrated values. Plots divide by bin width to show **flux density (per
keV)**, which is easier to interpret on irregular energy grids.

### Default test parameters

| Component | Parameters |
|-----------|------------|
| Gradus (lamppost) | spin=0.998, Eddington=0.1, inc=30°, h=3 r_g |
| xillver reflection | Γ=2.0, A_Fe=1.0, log ξ=2.0, n_e=17, incl=45° |

XSPEC adds `norm` and applies `*redshift` from `model.dat`; the Julia API uses
the nine physics parameters only.

## Figures

### 1. Line-profile kernel `L(g)`

Unit-area kernel in redshift space, `g = ν_obs / ν_em`.

![Line profile kernel](figures/01_line_profile_kernel.png)

**Expect:** a broadened profile peaked near `g ≈ 1`, falling to zero at the grid
edges, with ∫L(g) dg = 1.

### 2. Rest-frame reflection spectrum `R(E)`

Interpolated xillver spectrum on the table energy grid, plotted as flux density
(per keV) on a log–log scale.

![Reflection spectrum full](figures/02_reflection_spectrum.png)

### 3. Fe Kα region (3–9 keV)

Same reflection spectrum (flux per keV), linear scale, zoomed to the iron-line band.

![Reflection Fe region](figures/03_reflection_fe_region.png)

**Expect:** continuum plus Fe emission complex; exact shape depends on table
parameters.

### 4. Convolution with a narrow kernel

Reflection before and after convolution with a near-delta kernel (`g ≈ 1`). The
two curves should overlap closely in shape (peak location preserved).

![Narrow-kernel convolution](figures/04_convolution_narrow_kernel.png)

### 5. Convolution with the Gradus kernel

Reflection before and after blurring with the relativistic lamppost line profile.

![Gradus-kernel convolution](figures/05_convolution_gradus_kernel.png)

**Expect:** Fe features broadened and shifted — red wing extended, peak smeared.
This is the main science check before the full XSPEC model is wired up.

### 6. Flux-conserving rebin

Convolved spectrum on the table grid vs the same spectrum rebinned to a coarser
grid. Plotted as flux per keV; integrated flux per bin is conserved by the rebin.

![Rebin example](figures/06_rebin_example.png)

## Pipeline (current status)

```text
xillver FITS  →  R(E)              [table_model.jl]
Gradus params  →  L(g)              [line_profile.jl]
R, L           →  F = M * R         [convolution.jl]
F              →  XSPEC energies    [spectrum.jl: evaluate_spectrum]
```

`evaluate_spectrum` / `evaluate_spectrum_interpolated` orchestrate the full
pipeline in Julia. The XSPEC entry point (`gradusxspec`) calls the interpolated
version with nine physics parameters; `*redshift` and `norm` are handled by XSPEC.
