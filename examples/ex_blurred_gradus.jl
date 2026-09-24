#!/usr/bin/env julia
# Example 3: blurred reflection spectra for the Gradus corona models.
#
# Evaluates the full XSPEC pipeline (xillver table -> L(g) kernel -> blur ->
# rebin) via `evaluate_spectrum_interpolated`, the same entry point the XSPEC
# models call. Parameter values are chosen on the fitting grids so each scan
# point is a single grid corner: kernels are cached on disk
# (GRADUSXSPEC_KERNEL_CACHE), making re-runs fast.
#
# Usage (from repository root):
#   julia --project=. examples/ex_blurred_gradus.jl [path/to/xillverD-5.fits]
#
# First run computes every kernel with Gradus and can take tens of minutes
# (the gradus_disc_thin scan is the slowest: it stacks ring emissivities out
# to the corona outer radius). Subsequent runs use the kernel cache.

include(joinpath(@__DIR__, "common.jl"))

using GradusXSPEC

const REFL = (2.0, 1.0, 2.0, 17.0, 45.0)  # Γ, A_Fe, logξ, log n_e, incl
# The default blur working grid is tuned for fits above GRADUSXSPEC_BLUR_EMIN
# (2 keV); below that the grid coarsens quickly, so plot the 2-10 keV band.
const EDGES = collect(log_energy_edges(2.0, 10.0, 300))
const E_MID = 0.5 .* (EDGES[1:end-1] .+ EDGES[2:end])

function blurred_flux(model_name::AbstractString, gradus_params, table_path)
    rt = get_model_runtime(model_name)
    flux = evaluate_spectrum_interpolated(
        rt,
        EDGES,
        (gradus_params..., REFL...);
        table_path = table_path,
    )
    assert_finite_nonneg("$model_name$gradus_params", flux)
    return flux_per_keV(flux, EDGES[1:end-1], EDGES[2:end])
end

function scan_plot(model_name, curves, table_path; title, filename)
    p = plot(;
        xscale = :log10,
        yscale = :log10,
        xlabel = "Energy (keV)",
        ylabel = "Flux density (per keV)",
        title = title,
        legend = :bottomleft,
    )
    for (label, gradus_params) in curves
        println("  $model_name $label ...")
        F = blurred_flux(model_name, gradus_params, table_path)
        plot!(p, E_MID, positive_floor(F); label = label, linewidth = 1.5)
    end
    save_plot(p, filename)
end

function main()
    table_path = require_table(table_path_from_args())

    # Lamppost height: lower source -> stronger central illumination -> broader
    # red wing.
    println("Lamppost height scan (gradus_lamp_thin):")
    scan_plot(
        "gradus_lamp_thin",
        [("h = $h r_g", (0.998, 30.0, h)) for h in (2.0, 5.0, 10.0, 20.0)],
        table_path;
        title = "gradus_lamp_thin: lamppost height (a = 0.998, i = 30°)",
        filename = "ex3_blur_lamp_height.png",
    )

    # Spin: smaller ISCO at high spin extends the red wing.
    println("Spin scan (gradus_lamp_thin):")
    scan_plot(
        "gradus_lamp_thin",
        [("a = $a", (a, 30.0, 3.0)) for a in (0.1, 0.5, 0.9, 0.998)],
        table_path;
        title = "gradus_lamp_thin: spin (h = 3 r_g, i = 30°)",
        filename = "ex3_blur_lamp_spin.png",
    )

    # Inclination: Doppler shifts move the blue edge.
    println("Inclination scan (gradus_lamp_thin):")
    scan_plot(
        "gradus_lamp_thin",
        [("i = $(Int(inc))°", (0.998, inc, 3.0)) for inc in (15.0, 30.0, 45.0, 60.0)],
        table_path;
        title = "gradus_lamp_thin: inclination (a = 0.998, h = 3 r_g)",
        filename = "ex3_blur_lamp_incl.png",
    )

    # Thin disc vs Shakura–Sunyaev at the same lamppost geometry.
    println("Disc comparison (gradus_lamp_thin vs gradus_lamp_ss):")
    p = plot(;
        xscale = :log10,
        yscale = :log10,
        xlabel = "Energy (keV)",
        ylabel = "Flux density (per keV)",
        title = "Thin disc vs Shakura–Sunyaev (a = 0.998, i = 30°, h = 3 r_g)",
        legend = :bottomleft,
    )
    F_thin = blurred_flux("gradus_lamp_thin", (0.998, 30.0, 3.0), table_path)
    F_ss = blurred_flux("gradus_lamp_ss", (0.998, 0.1, 30.0, 3.0), table_path)
    plot!(p, E_MID, positive_floor(F_thin); label = "thin disc", linewidth = 1.5)
    plot!(p, E_MID, positive_floor(F_ss); label = "Shakura–Sunyaev (Ṁ_Edd = 0.1)", linewidth = 1.5)
    save_plot(p, "ex3_blur_disc_comparison.png")

    # Ring corona radius.
    println("Ring radius scan (gradus_ring_thin):")
    scan_plot(
        "gradus_ring_thin",
        [("r = $r r_g", (0.998, 30.0, r, 5.0)) for r in (2.0, 5.0, 10.0, 20.0)],
        table_path;
        title = "gradus_ring_thin: ring radius (a = 0.998, i = 30°, h = 5 r_g)",
        filename = "ex3_blur_ring_radius.png",
    )

    # Filled disc corona outer radius (slowest scan: R rings are stacked, but
    # rings are shared between R values through the ring-emissivity cache).
    println("Disc-corona outer radius scan (gradus_disc_thin):")
    scan_plot(
        "gradus_disc_thin",
        [("R = $R r_g", (0.998, 30.0, R, 5.0)) for R in (5.0, 10.0, 20.0)],
        table_path;
        title = "gradus_disc_thin: corona outer radius (a = 0.998, i = 30°, h = 5 r_g)",
        filename = "ex3_blur_disc_radius.png",
    )

    println("Done: 6 figures in ", OUT_DIR)
end

main()
