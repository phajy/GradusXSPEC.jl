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
| `kerrz_lamp_thin` | Lamppost ([kerrz](https://git.sr.ht/~fjebaker/kerrz) CLI) | Thin disc | spin, inc, h |
| `kerrz_ring_thin` | Co-rotating ring (kerrz CLI) | Thin disc | spin, inc, r, h |
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
  Implemented as a nested stack of ring coronae on a fixed radial mesh from
  1 r_g to the top of the `r` fitting grid, with the same spacing as that grid
  (all rings with radius ≤ `r`) and `r·Δr` weighting (uniform surface
  brightness). Because the mesh is independent of `r`, cached ring emissivity
  profiles are reused across outer radii.

### Kerrz CLI models

- **`kerrz_lamp_thin`** / **`kerrz_ring_thin`** — same XSPEC parameters as the
  Gradus thin lamppost / ring models, but L(g) is computed by shelling out to
  the [kerrz](https://git.sr.ht/~fjebaker/kerrz) CLI (`emissivity` → FITS →
  `lineprof`). Point `GRADUSXSPEC_KERRZ` at the binary if it is not at the
  default path. Useful for side-by-side comparison with `gradus_*_thin`.

### Diagnostic model

- **`test_gauss`** — convolves the reflection table with a Gaussian in
  `g = E_obs/E_em`. For very narrow `Sigma` the kernel is treated as an identity
  (no blur), useful for checking the table interpolation path without ray
  tracing.

## Parameter limits

These hard limits are enforced in `model.dat` to avoid known Gradus / kerrz
failure modes. They may be relaxed when fixed upstream.

| Parameter | Limit | Reason |
|-----------|-------|--------|
| `inc` | ≤ 65° | Gradus transfer-function failures at high inclination for some `(spin, h)` |
| `h` (ring/disc) | ≥ 2.5 r_g | Gradus `DomainError` in ring emissivity for low corona heights |
| `spin` | ≥ 0.1 | kerrz ring at spin = 0: empty disc emissivity → `lineprof` NaNError (see `kerrz_bugs.md`) |
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

## Blur working grid and soft X-ray band

Relativistic blurring is applied on a hybrid **working energy grid** (see
`src/blur_grid.jl`), not directly on the xillver table bins. By default the
**core band** runs from `GRADUSXSPEC_BLUR_EMIN` = **2 keV** to
`GRADUSXSPEC_BLUR_EMAX` = 150 keV with bin widths
`ΔE = max(0.1 keV, 1% × E)`. Below 2 keV the grid adds only a single extremely
coarse pad bin down to the table edge (and similarly above the core toward
high energies). That design targets Fe Kα and typical 2–10 keV fits without
paying the cost of a fine grid over the full xillver range.

**Consequence:** blurred model output is **zero below about 1.3 keV** with
default settings. Observed energies in that band fall into the coarse pad,
where the convolution cannot resolve reflected flux from the table; only once
the working grid reaches the fine core (near 2 keV) does structure appear.
This is expected, not a bug in the reflection table.

**Workarounds today:**

- Lower `GRADUSXSPEC_BLUR_EMIN` before starting XSPEC if you need blurred
  flux at softer energies (rebuild/restart not required — read at evaluation
  time).
- Set `GRADUSXSPEC_BLUR_NATIVE=1` to blur on the native xillver table grid
  instead (slower, but full table coverage).

**Possible future default:** extend the core band downward to about **0.1 keV**
with moderate resolution (similar `ΔE` rules as the current 2–150 keV core),
so soft-band and broad-band fits work out of the box without env overrides.
That would increase blur-grid size and memory use slightly; the trade-off has
not been implemented yet.

See also the `GRADUSXSPEC_BLUR_*` variables in [Using in XSPEC](xspec.md).

## Physics references

Implementation follows Gradus corona and disc types; see the
[Gradus line profiles documentation](https://astro-group.codeberg.page/Gradus.jl/dev/lineprofiles/)
and the [Gradus.jl source](https://codeberg.org/astro-group/Gradus.jl).

!!! note "Documentation status"
    Parameter descriptions, example spectra, and validation plots for the ring
    and disc corona models will be expanded after Linux build verification.
