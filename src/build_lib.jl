using PackageCompiler
include("model_definition.jl")

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

# Create the model.dat file to go along with the model.
write("model.dat", model_dat_text(; include_reflection = true))
