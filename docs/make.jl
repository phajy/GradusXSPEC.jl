using Documenter
using GradusXSPEC

makedocs(
    modules = [GradusXSPEC],
    sitename = "GradusXSPEC.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = [
        "Home" => "index.md",
        "Building" => "build.md",
        "Models" => "models.md",
        "Validation" => "validation.md",
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
)
