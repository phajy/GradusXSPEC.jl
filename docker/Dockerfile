# syntax=docker/dockerfile:1
#
# Linux build + smoke test for GradusXSPEC, layered on the official HEASoft
# container image.
#
# Prerequisite: build the HEASoft base image first, following
#   https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/docker.html
# That produces a local image tagged e.g. `heasoft:v6.36`. There is no official
# image on Docker Hub, so this must exist locally before building this one.
#
# Build (from the repository root, so the whole repo is the build context):
#   docker build -f docker/Dockerfile \
#     --build-arg HEASOFT_TAG=v6.36 \
#     --build-arg JULIA_VERSION=1.12.6 \
#     -t gradusxspec-linux .
#
# Run the smoke test (loads the package in XSPEC and evaluates a model):
#   docker run --rm gradusxspec-linux
#
# The reflection table (xillverD-5.fits) is gitignored and may not be in the
# build context. If it is absent, mount it at run time:
#   docker run --rm -v /path/to/xillverD-5.fits:/work/GradusXSPEC.jl/xillverD-5.fits \
#     gradusxspec-linux

ARG HEASOFT_TAG=v6.36
FROM heasoft:${HEASOFT_TAG}

ARG JULIA_VERSION=1.12.6

USER root

# Install a pinned Julia matching Project.toml. The HEASoft image is Ubuntu-based.
RUN set -eux; \
    arch="$(uname -m)"; \
    case "$arch" in \
      x86_64)  jarch=x64;      jsuffix=x86_64 ;; \
      aarch64) jarch=aarch64;  jsuffix=aarch64 ;; \
      *) echo "unsupported architecture: $arch" >&2; exit 1 ;; \
    esac; \
    minor="${JULIA_VERSION%.*}"; \
    url="https://julialang-s3.julialang.org/bin/linux/${jarch}/${minor}/julia-${JULIA_VERSION}-linux-${jsuffix}.tar.gz"; \
    curl -fsSL "$url" -o /tmp/julia.tar.gz; \
    tar -xzf /tmp/julia.tar.gz -C /opt; \
    ln -s "/opt/julia-${JULIA_VERSION}/bin/julia" /usr/local/bin/julia; \
    rm /tmp/julia.tar.gz; \
    julia --version

WORKDIR /work/GradusXSPEC.jl
COPY . .

RUN chmod +x docker/with-headas.sh *.sh

# Build the Julia shared library, then the XSPEC local-model package. HEADAS is
# sourced explicitly because Docker RUN steps do not inherit the interactive
# HEASoft environment.
RUN docker/with-headas.sh bash -c './build-julia.sh && ./build-xspec.sh'

# Default: run the headless smoke test.
CMD ["docker/with-headas.sh", "xspec", "-", "docker/smoke-test.xcm"]
