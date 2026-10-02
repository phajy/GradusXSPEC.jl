# Status — GradusXSPEC.jl

**One-liner:** Call Gradus.jl (and kerrz CLI) relativistic reflection models from XSPEC as local models.
**Collaborators:** Andy Young (maintainer); Gradus / kerrz upstream (Fergus Baker / Bristol astro group)
**Last updated:** 2026-10-02

## Where we are
- On `main` (synced with GitHub). Package version **0.2.0**. Six XSPEC models ship: `gradus_lamp_ss`, `gradus_lamp_thin`, `gradus_ring_thin`, `gradus_disc_thin`, `kerrz_lamp_thin`, `kerrz_ring_thin`.
- Performance path is in place: fixed disc-corona radial mesh + L1 ring cache, unlimited on-disk `L(g)` kernels, Float32 blur-matrix LRU (~16 GiB default), coarse 2–150 keV blur grid, **log-energy FFT convolution as default**.
- Build/docs: `./build-julia.sh` + `./build-xspec.sh`, Documenter manual, example gallery, validation scripts. Docs CI deploys to GitHub Pages ([stable](https://phajy.github.io/GradusXSPEC.jl/stable/) / [dev](https://phajy.github.io/GradusXSPEC.jl/dev/)); AstroRegistry added so Gradus resolves in CI/clones.
- Known upstream failures are worked around with parameter floors/ceilings (`spin ≥ 0.1`, `inc ≤ 65°`, ring/disc `h ≥ 2.5`) and a clamp on slightly negative ring-branch ε. Details live in untracked `gradus_bugs.md` / `kerrz_bugs.md` plus `reproduce_*.jl` scripts.
- Local kerrz work (hangs / race) was investigated; for now GradusXSPEC still shells out to the CLI rather than loading kerrz as a library.

## Next steps
- [ ] Merge or re-apply `fix/docs-ci-general-registry` (ensure General registry before AstroRegistry) so docs CI on `main` stops failing on unregistered deps (e.g. MutableArithmetics).
- [ ] File / chase upstream Gradus fixes: high-`inc` transfer functions, RingCorona sub-horizon `sqrt`, `DiscCorona` API + profile stacking, `log2` on ε≤0 (`gradus_bugs.md`).
- [ ] File / chase kerrz fixes: empty ring emissivity at `a = 0`; re-test lamppost hang after newer builds (`kerrz_bugs.md`). Longer-term: kerrz as a dylib instead of CLI (deferred until API settles).
- [ ] Optional model feature: expose corona photon index Γ as an XSPEC parameter (default 2, Δ=0.25, ~1–3) so it can be tied to `Refl_Gamma` (currently hardcoded for kerrz; Gradus path similarly fixed).
- [ ] Soft-band blur: consider extending the default core grid down to ~0.1 keV (documented as a future default; today output is ~zero below ~1.3 keV).
- [ ] Optional: `kerrz_disc_thin` once disc stacking / lineprof path is solid in kerrz (lamp + ring only for now).
- [ ] Expand ring/disc parameter docs and validation plots (noted in Documenter as pending after Linux build verification).
- [ ] Decide whether to track `gradus_bugs.md` / `kerrz_bugs.md` / `reproduce_*.jl`, or keep them personal/untracked.

## Open questions / blockers
- Docs CI on `main` after PR #4: General-registry gap may still break deploy until the fix branch is merged.
- Upstream Gradus `DiscCorona` remains unusable; GradusXSPEC keeps its own fixed-mesh stack.
- Fitting can still hit bad **grid corners** (exact continuous params OK) for Gradus transfer/emissivity and kerrz hangs — limits reduce but do not eliminate this.
- Soft X-ray fits need an explicit `GRADUSXSPEC_BLUR_EMIN` (or `BLUR_NATIVE`) until the default band is widened.

## Last session
- Date: 2026-09-24 (docs/CI); STATUS template filled 2026-10-02
- Did: Merged docs + examples + Documenter CI (PR #4); added AstroRegistry for CI/clones; earlier disc-mesh, FFT blur, kerrz backends, and TBB malloc stub are on `main`.
- Left off at: Docs CI registry fix branch exists but may not be on `main`; STATUS.md newly created; local untracked bug notes / repro scripts and a dirty `Manifest.toml` not committed.
