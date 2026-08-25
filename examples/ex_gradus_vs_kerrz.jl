#!/usr/bin/env julia
# Example 5: Gradus vs kerrz line profiles and blurred spectra.
#
# The kerrz_* models are parallel implementations of the lamppost and ring
# models that obtain L(g) from the external kerrz CLI instead of Gradus.jl.
# Two entirely independent ray tracers agreeing on the same physics is the
# strongest end-to-end check available, so this script overlays:
#   - unit-area L(g) kernels (lamppost and ring geometries),
#   - fully blurred reflection spectra through the XSPEC evaluation path.
#
# Requires the kerrz binary (default ~/GitHub/kerrz/zig-out/bin/kerrz, or set
# GRADUSXSPEC_KERRZ). Some difference between the codes is expected: they use
# different emissivity photon statistics, outer radii, and integration grids.
# The script fails only on gross disagreement (relative L2 > 0.5) and prints
# the actual values for the documentation.
#
# Usage (from repository root, ~minutes):
#   julia --project=. examples/ex_gradus_vs_kerrz.jl [path/to/xillverD-5.fits]

include(joinpath(@__DIR__, "common.jl"))

using LinearAlgebra
using GradusXSPEC

const REFL = (2.0, 1.0, 2.0, 17.0, 45.0)
# Match ex_blurred_gradus.jl: the default blur grid targets >= 2 keV.
const EDGES = collect(log_energy_edges(2.0, 10.0, 300))
const E_MID = 0.5 .* (EDGES[1:end-1] .+ EDGES[2:end])

rel_l2(a, b) = norm(a - b) / max(norm(b), eps())

function kernel_comparison(params, gradus_variant, kerrz_variant; title, filename)
    println("  Gradus $gradus_variant $params ...")
    g, L_gradus = line_profile_kernel(params, gradus_variant, :thin)
    println("  kerrz $kerrz_variant $params ...")
    _, L_kerrz = line_profile_kernel(params, kerrz_variant, :thin)
    rel = rel_l2(L_kerrz, L_gradus)
    println("  kernel relative L2: ", rel)
    rel <= 0.5 || error("Gradus and kerrz kernels grossly disagree (L2 = $rel) for $params")

    p_top = plot(
        g,
        L_gradus;
        ylabel = "L(g) (unit area)",
        label = "Gradus",
        linewidth = 2,
        title = title * "  (relative L2 = $(round(rel, sigdigits = 2)))",
    )
    plot!(p_top, g, L_kerrz; label = "kerrz", linewidth = 2, linestyle = :dash)
    p_bottom = plot(
        g,
        L_kerrz .- L_gradus;
        xlabel = "g = E_obs / E_em",
        ylabel = "kerrz − Gradus",
        label = "",
        linewidth = 1,
    )
    hline!(p_bottom, [0.0]; color = :gray, linestyle = :dash, label = "")
    p = plot(p_top, p_bottom; layout = (2, 1), link = :x, size = (700, 600))
    save_plot(p, filename)
    return rel
end

function blurred_comparison(gradus_model, kerrz_model, gradus_params, table_path; title, filename)
    println("  $gradus_model $gradus_params ...")
    rt_g = get_model_runtime(gradus_model)
    F_gradus = evaluate_spectrum_interpolated(
        rt_g, EDGES, (gradus_params..., REFL...); table_path = table_path)
    println("  $kerrz_model $gradus_params ...")
    rt_k = get_model_runtime(kerrz_model)
    F_kerrz = evaluate_spectrum_interpolated(
        rt_k, EDGES, (gradus_params..., REFL...); table_path = table_path)
    assert_finite_nonneg(gradus_model, F_gradus)
    assert_finite_nonneg(kerrz_model, F_kerrz)
    rel = rel_l2(F_kerrz, F_gradus)
    println("  blurred-spectrum relative L2: ", rel)
    rel <= 0.5 || error("Blurred spectra grossly disagree (L2 = $rel)")

    lo, hi = EDGES[1:end-1], EDGES[2:end]
    Fg = flux_per_keV(F_gradus, lo, hi)
    Fk = flux_per_keV(F_kerrz, lo, hi)
    denom = max.(Fg, 1e-12 * maximum(Fg))

    p_top = plot(
        E_MID,
        positive_floor(Fg);
        xscale = :log10,
        yscale = :log10,
        ylabel = "Flux density (per keV)",
        label = gradus_model,
        linewidth = 2,
        title = title * "  (relative L2 = $(round(rel, sigdigits = 2)))",
    )
    plot!(p_top, E_MID, positive_floor(Fk);
        label = kerrz_model, linewidth = 2, linestyle = :dash)
    p_bottom = plot(
        E_MID,
        (Fk .- Fg) ./ denom;
        xscale = :log10,
        xlabel = "Energy (keV)",
        ylabel = "(kerrz − Gradus) / Gradus",
        label = "",
        linewidth = 1,
    )
    hline!(p_bottom, [0.0]; color = :gray, linestyle = :dash, label = "")
    p = plot(p_top, p_bottom; layout = (2, 1), link = :x, size = (700, 600))
    save_plot(p, filename)
    return rel
end

function main()
    bin = try
        kerrz_binary_path()
    catch e
        error(
            "kerrz binary not found ($(sprint(showerror, e))).\n" *
            "Build kerrz (zig build -Doptimize=ReleaseFast in the kerrz repository)\n" *
            "or set GRADUSXSPEC_KERRZ to the executable.",
        )
    end
    println("Using kerrz binary: ", bin)
    table_path = require_table(table_path_from_args())

    println("Lamppost kernel comparison (a = 0.998, i = 30°, h = 5 r_g):")
    kernel_comparison(
        (0.998, 30.0, 5.0), :lamppost, :kerrz_lamppost;
        title = "Lamppost L(g): Gradus vs kerrz",
        filename = "ex5_kernel_lamppost.png",
    )

    println("Ring kernel comparison (a = 0.998, i = 30°, r = 5 r_g, h = 5 r_g):")
    kernel_comparison(
        (0.998, 30.0, 5.0, 5.0), :ring, :kerrz_ring;
        title = "Ring L(g): Gradus vs kerrz",
        filename = "ex5_kernel_ring.png",
    )

    println("Blurred spectrum comparison (lamppost):")
    blurred_comparison(
        "gradus_lamp_thin", "kerrz_lamp_thin", (0.998, 30.0, 5.0), table_path;
        title = "Blurred reflection: Gradus vs kerrz lamppost",
        filename = "ex5_blurred_lamppost.png",
    )

    println("Done: 3 figures in ", OUT_DIR)
end

main()
