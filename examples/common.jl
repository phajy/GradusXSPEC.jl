# Shared helpers for the example-gallery scripts (examples/ex_*.jl).
#
# Figures are written to docs/src/figures/ so the Documenter manual can embed
# them. Each script also performs cheap numeric sanity checks so a broken
# pipeline fails loudly instead of producing silently wrong plots.

try
    using Plots
catch
    error(
        "Plots.jl is required. Run once:\n" *
        "  julia --project=. -e 'using Pkg; Pkg.add(\"Plots\")'",
    )
end

const OUT_DIR = joinpath(@__DIR__, "..", "docs", "src", "figures")
const DEFAULT_TABLE = "xillverD-5.fits"

table_path_from_args() = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_TABLE

function require_table(path::AbstractString)
    isfile(path) || error(
        "table model not found: $path\n" *
        "Download xillverD-5.fits (see docs) or pass a path as the first argument.",
    )
    return path
end

function save_plot(p, filename::String)
    mkpath(OUT_DIR)
    path = joinpath(OUT_DIR, filename)
    savefig(p, path)
    println("Wrote ", path)
end

energy_midpoints(lo, hi) = (Float64.(lo) .+ Float64.(hi)) ./ 2
bin_widths(lo, hi) = Float64.(hi) .- Float64.(lo)

"""XSPEC-style bin-integrated flux -> flux density (per keV), for plotting."""
flux_per_keV(flux, lo, hi) = Float64.(flux) ./ bin_widths(lo, hi)

"""Trapezoidal integral of samples `y` on grid `x`."""
function integrate_grid(x, y)
    area = 0.0
    for i in 1:(length(x) - 1)
        area += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return area
end

"""Fail loudly if `v` contains non-finite or (beyond tolerance) negative values."""
function assert_finite_nonneg(name::AbstractString, v; neg_tol::Float64 = 0.0)
    all(isfinite, v) || error("$name contains non-finite values")
    lo = minimum(v)
    lo >= -neg_tol || error("$name contains negative values (min=$lo)")
    return nothing
end

"""Logarithmically spaced energy bin edges (keV)."""
log_energy_edges(e_min, e_max, n_bins) =
    exp.(range(log(Float64(e_min)), log(Float64(e_max)), length = n_bins + 1))

"""Floor values for log-scale plotting without distorting real features."""
positive_floor(v; floor_frac = 1e-10) = max.(v, floor_frac * maximum(v))
