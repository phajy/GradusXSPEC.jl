# Visual validation (development)

This content now lives in the Documenter manual:

- [Validation](src/validation.md) — scripts, expected figures, pipeline
- [Building](src/build.md) — build and platform notes

Build and open the manual locally:

```sh
./build-docs.sh
open docs/build/index.html   # macOS
```

The scripts and figure paths are unchanged; `plot_validation.jl` still writes to
`docs/figures/`.
