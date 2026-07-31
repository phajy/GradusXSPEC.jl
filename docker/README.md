# Testing the Linux build with Docker

The macOS and Linux builds share the same pipeline (`build-julia.sh` →
`initpackage` → Makefile patch → `hmake`); only the gcc@14 Fortran-path
workaround is macOS-specific and is skipped elsewhere.

**If you already have HEASOFT installed on Linux, use the native build path
instead** — see the [Building](https://github.com/phajy/GradusXSPEC.jl/blob/main/docs/src/build.md#linux-native-heasoft)
section in the manual (`./check-env.sh`, then `./build-julia.sh` and
`./build-xspec.sh`). Docker is optional: it provides a reproducible clean-room
check when you do not want to rely on a particular host install.

## Prerequisites

1. **Docker** installed and running.
2. **A local HEASoft base image.** HEASARC does not publish an image on Docker
   Hub; you build one from their source tarball by following
   [HEASoft & Docker](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/docker.html).
   That produces a local image tagged like `heasoft:v6.36`.
3. **The reflection table** `xillverD-5.fits`. It is gitignored, so either place
   it in the repository root before building (it will be copied into the image)
   or mount it at run time (see below).

## Build

From the **repository root** (the whole repo is the build context):

```sh
docker build -f docker/Dockerfile \
  --build-arg HEASOFT_TAG=v6.36 \
  --build-arg JULIA_VERSION=1.12.6 \
  -t gradusxspec-linux .
```

This installs a pinned Julia, compiles the Julia shared library, and builds the
XSPEC local-model package inside the container. A build failure here is exactly
the Linux build regression we want to catch.

## Run the smoke test

```sh
docker run --rm gradusxspec-linux
```

This loads `gradusxspec` in XSPEC and evaluates the fast `test_gauss` diagnostic
model over a dummy response, printing a flux value. It confirms the XSPEC ↔
Julia bridge links on Linux. The physical (ray-traced) models use the same
interface; `docker/smoke-test.xcm` shows how to exercise one, but each
evaluation takes a couple of minutes.

If the table is not baked into the image, mount it:

```sh
docker run --rm \
  -v /path/to/xillverD-5.fits:/work/GradusXSPEC.jl/xillverD-5.fits \
  gradusxspec-linux
```

## Interactive shell

To poke around inside the container (HEADAS pre-sourced):

```sh
docker run --rm -it gradusxspec-linux docker/with-headas.sh bash
```

## GitHub Actions

A full Docker-based CI job is feasible but heavy (~12 GB image, long
PackageCompiler build). It is best run on `workflow_dispatch` or release tags
rather than every commit. See the [Continuous integration](https://github.com/phajy/GradusXSPEC.jl/blob/main/docs/src/build.md#continuous-integration)
section in the manual for a tiered CI strategy.

## Notes

- A segmentation fault on `quit`/`exit` after evaluations is a known, harmless
  interaction between the embedded Julia runtime and XSPEC teardown; it happens
  after results are produced.
- The `HEASOFT_TAG` and `JULIA_VERSION` build args let you match your local
  HEASoft image tag and the Julia version pinned in `Project.toml`.
