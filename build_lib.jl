using PackageCompiler

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

# Create the model.dat file to go along with the model
open("model.dat", "w") do f
    write(f, "gradus 1 0. 1.e20 c_gradusjulia add 0 0\n")
    write(f, "inc degrees 30.0 1.0 1.0 89.0 89.0 0.1\n")
end
