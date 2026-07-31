using PackageCompiler
include("model_definition.jl")
include(joinpath(dirname(@__DIR__), "scripts", "stub_tbbmalloc_proxy.jl"))

# This generates the library in a "build" folder
# It creates libGradusXSPEC.so (Linux) or .dll (Windows) and a header file
create_library(
    ".",                      # Path to the package we just made
    "build",                  # Output directory
    lib_name="libGradusXSPEC", 
    force=true,
    # header_files = ["src/GradusXSPEC_api.h"], # Optional: if you want custom headers, otherwise PC generates one
    include_transitive_dependencies=true,
    include_lazy_artifacts=true
)

# macOS: neutralise oneTBB's malloc-zone hijack, which segfaults XSPEC.
# See scripts/stub_tbbmalloc_proxy.jl for the full rationale.
stub_tbbmalloc_proxy!("build")

# Create model.dat for all XSPEC models in this package.
write("model.dat", model_dat_text())
