# Example gallery

Worked examples exercising each stage of the pipeline, from the un-blurred
reflection table to full Gradus-vs-kerrz cross-checks. Every figure on this
page is produced by a script in
[`examples/`](https://github.com/phajy/GradusXSPEC.jl/tree/main/examples);
regenerate them all with:

```sh
julia --project=. examples/ex_reflection.jl        # seconds
julia --project=. examples/ex_line_profiles.jl     # minutes
julia --project=. examples/ex_convolution.jl       # ~1-2 min
julia --project=. examples/ex_blurred_gradus.jl    # tens of minutes cold, fast cached
julia --project=. examples/ex_gradus_vs_kerrz.jl   # minutes; needs kerrz binary
```

Each script asserts basic invariants (finite non-negative fluxes, unit kernel
areas, agreement tolerances) so a broken pipeline fails loudly rather than
producing quietly wrong figures. Figures are written to `docs/src/figures/`
(gitignored); rebuild this manual with `./build-docs.sh` to embed them.

## 1. Un-blurred reflection spectra

`examples/ex_reflection.jl` interpolates the xillver table at a range of
values for each parameter, holding the others at Γ = 2, A_Fe = 1, log ξ = 2,
log n_e = 17, incl = 45°. These are the rest-frame spectra every model starts
from — no relativistic effects yet.

![Photon index scan](figures/ex1_refl_gamma.png)

Steeper Γ softens the continuum; the fluorescent lines ride on top.

![Iron abundance scan](figures/ex1_refl_afe.png)
![Iron abundance, Fe Kα zoom](figures/ex1_refl_afe_zoom.png)

Iron abundance mostly rescales the Fe Kα complex at 6.4–7 keV and the iron
K edge, which the zoom panel shows linearly.

![Ionisation scan](figures/ex1_refl_logxi.png)

Ionisation is the biggest lever: near-neutral discs (log ξ ≈ 0) show strong
narrow features and a hard reflection hump, while highly ionised discs
approach a featureless mirror.

![Density scan](figures/ex1_refl_dens.png)

Higher densities raise the soft-excess continuum below ~2 keV.

![Table inclination scan](figures/ex1_refl_incl.png)

The xillver inclination changes limb brightening of the reflected continuum;
note this is the *table* viewing angle, distinct from the ray-traced disc
inclination in the blurring.

## 2. Disc lines with power-law emissivity

`examples/ex_line_profiles.jl` computes the classic diskline test case: a thin
disc whose emissivity is an analytic power law ε(r) ∝ r^(-q), traced with the
same Gradus `lineprofile` call the models use. There is no corona here — this
isolates the transfer-function integration so the profile shapes can be
checked against intuition independently of any emissivity calculation. (This
configuration is deliberately not exported to XSPEC; it exists to build
confidence in the convolution kernels.)

![Emissivity index scan](figures/ex2_line_q_scan.png)

Steeper emissivity (larger q) weights the innermost radii, extending the
gravitationally redshifted wing to lower g.

![Spin scan](figures/ex2_line_spin_scan.png)

Higher spin moves the ISCO inward — from 4.23 r_g at a = 0.5 to 1.24 r_g at
a = 0.998 — feeding the red wing.

![Inclination scan](figures/ex2_line_incl_scan.png)

Inclination controls the Doppler blue edge: nearly face-on discs give narrow
profiles near g ≈ 1, edge-on discs push the blue horn beyond g > 1.2.

## 3. Blurred reflection: the Gradus models

`examples/ex_blurred_gradus.jl` runs the full XSPEC pipeline
(`evaluate_spectrum_interpolated`, the same entry point the XSPEC wrappers
call) for each model. Scan values sit on the fitting grids, so each curve is
a single cached grid corner. Spectra are shown over 2–10 keV, where the
default blur grid resolves structure reliably; see
[Models — blur grid](models.md#blur-working-grid-and-soft-x-ray-band) for why
output is zero below ~1.3 keV and how that might change in a future release.

![Lamppost height scan](figures/ex3_blur_lamp_height.png)

A lower lamppost concentrates illumination on the inner disc: the Fe line
red wing broadens and the features smear together.

![Spin scan](figures/ex3_blur_lamp_spin.png)
![Inclination scan](figures/ex3_blur_lamp_incl.png)

The blurred spectra inherit the line-profile behaviour from example 2: spin
extends the red wing, inclination moves the blue edge of every feature.

![Thin vs Shakura–Sunyaev disc](figures/ex3_blur_disc_comparison.png)

`gradus_lamp_ss` traces a finite-thickness Shakura–Sunyaev disc; at
Ṁ/Ṁ_Edd = 0.1 the differences from the razor-thin disc are subtle changes in
the innermost illumination and shadowing.

![Ring radius scan](figures/ex3_blur_ring_radius.png)

A ring corona at growing radius r illuminates progressively larger disc radii,
narrowing the blurring. At large r and h the profile approaches the lamppost
limit.

![Disc-corona outer radius scan](figures/ex3_blur_disc_radius.png)

The filled disc corona stacks ring contributions out to R (uniform surface
brightness, mesh spacing 1 r_g); increasing R adds weakly blurred flux from
large radii on top of the strongly blurred core.

## 4. Convolution computed in different ways

`examples/ex_convolution.jl` blurs the same xillver spectrum with the same
Gradus lamppost kernel through both convolution implementations — the
bin-integrated matrix path and the default log-energy FFT path — and shows the
overlay, the relative difference, and timings. The two implementations share
no code, so sub-percent agreement checks both. The script errors if the
relative L2 difference exceeds 1%.

![Matrix vs FFT convolution](figures/ex4_convolution_methods.png)

The FFT path is the default (`GRADUSXSPEC_CONVOLVE=fft`) and is substantially
faster once warm; the matrix path remains available with
`GRADUSXSPEC_CONVOLVE=matrix`.

## 5. Gradus vs kerrz

`examples/ex_gradus_vs_kerrz.jl` compares the two ray-tracing backends —
Gradus.jl and the [kerrz](https://git.sr.ht/~fjebaker/kerrz) CLI — on the same
geometries. These are fully independent codes (different language, integrator,
emissivity method), so agreement is the strongest end-to-end validation
available. Exact agreement is not expected: the backends use different photon
statistics, outer radii, and integration grids; the residual panels quantify
the difference, and the script fails only on gross (>50% L2) disagreement.

![Lamppost kernel comparison](figures/ex5_kernel_lamppost.png)
![Ring kernel comparison](figures/ex5_kernel_ring.png)

Unit-area L(g) kernels for the lamppost (a = 0.998, i = 30°, h = 5 r_g) and
ring (r = h = 5 r_g) geometries.

![Blurred spectrum comparison](figures/ex5_blurred_lamppost.png)

The same comparison propagated through the full pipeline:
`gradus_lamp_thin` vs `kerrz_lamp_thin` at identical parameters. Differences
in the kernels are diluted by the convolution with the broad reflection
continuum.
