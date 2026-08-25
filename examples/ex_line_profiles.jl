#!/usr/bin/env julia
# Example 2: simple disc-line profiles with power-law emissivity.
#
# Computes relativistic line profiles L(g) with Gradus for a thin disc whose
# emissivity is a plain power law ε(r) ∝ r^(-q) — the classic diskline/kerrline
# setup. There is no corona model here: this isolates the transfer-function
# integration that all GradusXSPEC models rely on, so the shapes can be checked
# against intuition (and the literature) independently of any emissivity
# calculation. This model is intentionally not exported to XSPEC.
#
# Uses the same Gradus call as src/line_profile.jl (`lineprofile` with
# TransferFunctionMethod and maxrₑ = 400), only with the analytic emissivity in
# place of a corona-derived profile.
#
# Usage (from repository root, ~minutes):
#   julia --project=. examples/ex_line_profiles.jl

include(joinpath(@__DIR__, "common.jl"))

using Gradus

const G_GRID = collect(range(0.1, 1.5, length = 300))

"""
    powerlaw_line_profile(; a, incl_deg, q) -> (g, L)

Unit-area line profile for a thin disc with ε(r) ∝ r^(-q), spin `a`, observer
inclination `incl_deg`.
"""
function powerlaw_line_profile(; a::Float64, incl_deg::Float64, q::Float64)
    m = KerrMetric(M = 1.0, a = a)
    d = ThinDisc(0.0, Inf)
    x = SVector(0.0, 1000.0, deg2rad(incl_deg), 0.0)
    bins = copy(G_GRID)
    _, flux = lineprofile(
        bins,
        r -> r^(-q),
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        maxrₑ = 400.0,
    )
    flux = Vector{Float64}(flux)
    assert_finite_nonneg("L(a=$a, incl=$incl_deg, q=$q)", flux; neg_tol = 1e-8)
    clamp!(flux, 0.0, Inf)
    area = integrate_grid(bins, flux)
    area > 0 || error("line profile integral is zero (a=$a, incl=$incl_deg, q=$q)")
    return bins, flux ./ area
end

function scan_plot(curves; title::AbstractString, filename::AbstractString)
    p = plot(;
        xlabel = "g = E_obs / E_em",
        ylabel = "L(g) (unit area)",
        title = title,
        legend = :topleft,
    )
    for (label, kwargs) in curves
        println("  computing $label ...")
        g, L = powerlaw_line_profile(; kwargs...)
        plot!(p, g, L; label = label, linewidth = 2)
    end
    vline!(p, [1.0]; linestyle = :dash, color = :gray, label = "")
    save_plot(p, filename)
end

function main()
    println("Emissivity-index scan (spin = 0.998, incl = 30°):")
    scan_plot(
        [("q = $q", (a = 0.998, incl_deg = 30.0, q = q)) for q in (2.0, 3.0, 4.0)];
        title = "Disc line, ε ∝ r^-q (a = 0.998, i = 30°)",
        filename = "ex2_line_q_scan.png",
    )

    println("Spin scan (q = 3, incl = 30°):")
    scan_plot(
        [("a = $a", (a = a, incl_deg = 30.0, q = 3.0)) for a in (0.5, 0.9, 0.998)];
        title = "Disc line, spin dependence (q = 3, i = 30°)",
        filename = "ex2_line_spin_scan.png",
    )

    println("Inclination scan (q = 3, spin = 0.998):")
    scan_plot(
        [("i = $(Int(i))°", (a = 0.998, incl_deg = i, q = 3.0)) for i in (15.0, 30.0, 45.0, 60.0)];
        title = "Disc line, inclination dependence (q = 3, a = 0.998)",
        filename = "ex2_line_incl_scan.png",
    )

    println("Done: 3 figures in ", OUT_DIR)
end

main()
