#!/usr/bin/env julia
# Example 1: un-blurred xillver reflection spectra for different parameters.
#
# Interpolates the reflection table at a range of values for each xillver
# parameter (photon index, iron abundance, ionisation, density, inclination)
# with everything else held at the defaults. No relativistic blurring is
# applied, so these are the rest-frame spectra the gradus_*/kerrz_* models
# start from.
#
# Usage (from repository root, ~seconds):
#   julia --project=. examples/ex_reflection.jl [path/to/xillverD-5.fits]

include(joinpath(@__DIR__, "common.jl"))

using GradusXSPEC

const BASE = (Gamma = 2.0, A_Fe = 1.0, logXi = 2.0, Dens = 17.0, Incl = 45.0)

refl_tuple(nt) = (nt.Gamma, nt.A_Fe, nt.logXi, nt.Dens, nt.Incl)

function scan_plot(
    table,
    param::Symbol,
    values;
    label::AbstractString,
    filename::AbstractString,
    e_range = (0.3, 100.0),
    logy::Bool = true,
)
    E = energy_midpoints(table.energy_lo, table.energy_hi)
    idx = findall(e_range[1] .<= E .<= e_range[2])
    p = plot(;
        xscale = :log10,
        yscale = logy ? :log10 : :identity,
        xlabel = "Energy (keV)",
        ylabel = "Flux density (per keV)",
        title = "xillver reflection: varying $label",
        legend = :bottomleft,
    )
    for v in values
        nt = merge(BASE, NamedTuple{(param,)}((Float64(v),)))
        R = interpolate_table_spectrum(table, refl_tuple(nt))
        assert_finite_nonneg("R($param=$v)", R)
        density = flux_per_keV(R, table.energy_lo, table.energy_hi)
        y = logy ? positive_floor(density[idx]) : density[idx]
        plot!(p, E[idx], y; label = "$label = $v", linewidth = 1.5)
    end
    save_plot(p, filename)
end

function main()
    table_path = require_table(table_path_from_args())
    println("Loading table: ", table_path)
    table = load_xspec_table(table_path)
    println("Base parameters: ", BASE)

    scan_plot(table, :Gamma, [1.6, 2.0, 2.4, 3.0];
        label = "Γ", filename = "ex1_refl_gamma.png")
    scan_plot(table, :A_Fe, [0.5, 1.0, 3.0, 10.0];
        label = "A_Fe", filename = "ex1_refl_afe.png")
    # Fe Kα zoom for the abundance scan: the line strength is the clearest effect.
    scan_plot(table, :A_Fe, [0.5, 1.0, 3.0, 10.0];
        label = "A_Fe", filename = "ex1_refl_afe_zoom.png",
        e_range = (3.0, 9.0), logy = false)
    scan_plot(table, :logXi, [0.0, 1.0, 2.0, 3.0, 4.0];
        label = "log ξ", filename = "ex1_refl_logxi.png")
    scan_plot(table, :Dens, [15.0, 17.0, 19.0];
        label = "log n_e", filename = "ex1_refl_dens.png")
    scan_plot(table, :Incl, [20.0, 45.0, 80.0];
        label = "incl (deg)", filename = "ex1_refl_incl.png")

    println("Done: 6 figures in ", OUT_DIR)
end

main()
