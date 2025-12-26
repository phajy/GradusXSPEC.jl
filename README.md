# GradusXSPEC.jl

Use Gradus to evaluate models from XSPEC

## How to build the library and import it into XSPEC

1. Make the Julia library

```sh
julia --project=. build_lib.jl
```

This will build the shared library `build/lib/libGradusXSPEC.dylib` and create the `model.dat` file.

2. Create the XSPEC model

```sh
initpackage gradus model.dat .
patch Makefile < changes.patch
hmake
```

3. Test the model

```
xspec
lmod gradus .
```

Set the model parameters, evaluate over a dummy energy grid, and use `iplot` to visualize the results.

```
========================================================================
Model gradus<1> Source No.: 1   Active/Off
Model Model Component  Parameter  Unit     Value
 par  comp
   1    1   gradus     spin                0.900000     +/-  0.0
   2    1   gradus     Eddington           0.150000     +/-  0.0
   3    1   gradus     inc        degrees  40.0000      +/-  0.0
   4    1   gradus     h          r_g      6.00000      +/-  0.0
   5    1   gradus     norm                1.00000      +/-  0.0
________________________________________________________________________
```

![XSPEC plot output](figs/pgplot.png)
