using Documenter
using GradusXSPEC

makedocs(
    modules = [GradusXSPEC],
    sitename = "GradusXSPEC.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Building" => "build.md",
        "Using in XSPEC" => "xspec.md",
        "Models" => "models.md",
        "Examples" => "examples.md",
        "Validation" => "validation.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/phajy/GradusXSPEC.jl.git",
    push_preview = true,
)
