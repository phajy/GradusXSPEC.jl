# GradusXSPEC.jl

Use Gradus to evaluate models from XSPEC

## How to build the library and import it into XSPEC

1. Make the Julia library

`julia --project=. build_lib.jl` in the terminal will build the shared library and create the `model.dat` file.

2. Create the XSPEC model

```sh
initpackage gradus model.dat .
patch Makefile < changes.patch
hmake

3. Test the model

`xspec`
`lmod gradus .`
