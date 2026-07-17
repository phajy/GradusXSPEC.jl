# Models

Each model convolves an xillver reflection table (`xillverD-5.fits` by default)
with a relativistic line profile computed in Gradus. All production models share
the same five reflection parameters (`Refl_Gamma`, `Refl_A_Fe`, `Refl_logXi`,
`Refl_Dens`, `Refl_Incl`) plus XSPEC's `norm`.

## Model summary

| XSPEC name | Corona | Disc | Gradus parameters |
|------------|--------|------|-------------------|
| `gradus_lamp_ss` | Lamppost | Shakura–Sunyaev (thick) | spin, Eddington, inc, h |
| `gradus_lamp_thin` | Lamppost | Thin disc | spin, inc, h |
| `gradus_ring_thin` | Co-rotating ring | Thin disc | spin, inc, r, h |
| `gradus_disc_thin` | Filled disc corona | Thin disc | spin, inc, r (outer), h |
| `test_gauss` | Gaussian blur in `g` | — | Sigma (+ reflection params) |

### Lamppost models

- **`gradus_lamp_ss`** — standard lamppost over a Shakura–Sunyaev thick disc.
  Includes an `Eddington` parameter controlling the disc scale height.
- **`gradus_lamp_thin`** — lamppost over a Novikov–Thorne thin disc (no
  `Eddington` parameter).

### Ring and disc corona models

- **`gradus_ring_thin`** — emission from a single co-rotating ring corona at
  radius `r` and height `h` above a thin disc.
- **`gradus_disc_thin`** — filled disc corona of outer radius `r` at height `h`.
  Implemented as a stack of ring coronae with `r·Δr` weighting (uniform surface
  brightness).

### Diagnostic model

- **`test_gauss`** — convolves the reflection table with a Gaussian in
  `g = E_obs/E_em`. For very narrow `Sigma` the kernel is treated as an identity
  (no blur), useful for checking the table interpolation path without ray
  tracing.

## Parameter limits (Gradus)

These hard limits are enforced in `model.dat` to avoid known Gradus failure modes.
They may be relaxed upstream in Gradus.jl later.

| Parameter | Limit | Reason |
|-----------|-------|--------|
| `inc` | ≤ 65° | Transfer-function failures at high inclination for some `(spin, h)` |
| `h` (ring/disc) | ≥ 2.5 r_g | `DomainError` in ring emissivity for low corona heights |
| `spin` | ≤ 0.998 | Standard Kerr bound |

Reflection parameters follow the xillver table ranges configured in
`model_definition.jl`.

## Init string

The model init string starts with the reflection table path (default
`xillverD-5.fits`). Optional tokens:

| Token | Effect |
|-------|--------|
| `verbose` | Print per-evaluation diagnostics |
| `monitor` | Write fit monitor file (default path) |
| `monitor=/path/to/file` | Custom monitor file path |
| `monitor_interval=N` | Refresh monitor every N evaluations |

Example: `xillverD-5.fits verbose monitor`

## Physics references

Implementation follows Gradus corona and disc types; see the
[Gradus line profiles documentation](https://astro-group.codeberg.page/Gradus.jl/dev/lineprofiles/)
and the [Gradus.jl source](https://codeberg.org/astro-group/Gradus.jl).

!!! note "Documentation status"
    Parameter descriptions, example spectra, and validation plots for the ring
    and disc corona models will be expanded after Linux build verification.
