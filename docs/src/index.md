# GradusXSPEC.jl

GradusXSPEC connects [Gradus.jl](https://astro-group.codeberg.page/Gradus.jl/dev/)
relativistic line profiles to [XSPEC](https://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/)
local models. Each model convolves an xillver reflection table with a Gradus
corona geometry and accretion-disc setup.

## Available models

| XSPEC name | Corona | Disc |
|------------|--------|------|
| `gradus_lamp_ss` | Lamppost | Shakura–Sunyaev (thick) |
| `gradus_lamp_thin` | Lamppost | Thin disc |
| `gradus_ring_thin` | Ring | Thin disc |
| `gradus_disc_thin` | Filled disc | Thin disc |
| `test_gauss` | Gaussian blur (diagnostic) | — |

Load the package once in XSPEC with `lmod gradusxspec .`, then select any model
name in `model`.

## Quick start

1. [Build the Julia library and XSPEC package](build.md).
2. [Check parameter meanings and limits](models.md).
3. [Validate spectra and generate plots](validation.md) before or alongside
   XSPEC testing.

## References

- [Appendix C: Adding Models to XSPEC](https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSappendixLocal.html)
- [OGIP 92-009 table model format (PDF)](https://heasarc.gsfc.nasa.gov/FTP/caldb/docs/memos/ogip_92_009/ogip_92_009.pdf)
- [Gradus.jl codebase](https://codeberg.org/astro-group/Gradus.jl)
- [Gradus.jl documentation](https://astro-group.codeberg.page/Gradus.jl/dev/)
