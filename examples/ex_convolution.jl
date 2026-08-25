#!/usr/bin/env julia
# Example 4: the same convolution computed in different ways.
#
# GradusXSPEC can apply the relativistic blur two ways:
#   - matrix: bin-integrated convolution matrix (legacy path,
#     GRADUSXSPEC_CONVOLVE=matrix)
#   - fft:    log-energy FFT via the convolution theorem (default,
#     GRADUSXSPEC_CONVOLVE=fft)
# The two are independent implementations of the same operator, so their
# agreement is a strong check on both. This script blurs the same xillver
# spectrum with the same Gradus lamppost kernel through both paths, plots the
# overlay and the relative difference, and reports timings.
#
# See scripts/compare_convolution_methods.jl for the plain numeric version
# used in validation (no plotting dependencies).
#
# Usage (from repository root, ~1–2 min; most of it is the Gradus kernel):
#   julia --project=. examples/ex_convolution.jl [path/to/xillverD-5.fits]

include(joinpath(@__DIR__, "common.jl"))

using LinearAlgebra
using GradusXSPEC

const REFL = (2.0, 1.0, 2.0, 17.0, 45.0)
const LAMP_PARAMS = (0.998, 30.0, 3.0)

function main()
    table_path = require_table(table_path_from_args())
    println("Loading table: ", table_path)
    table = load_xspec_table(table_path)
    R_tab = interpolate_table_spectrum(table, REFL)

    # FFT's native setting is a uniform log-energy grid; compare there.
    log_lo, log_hi = log_energy_bin_edges(0.5, 20.0, 512)
    R = rebin_flux(R_tab, table.energy_lo, table.energy_hi, log_lo, log_hi)

    println("Computing Gradus lamppost kernel ", LAMP_PARAMS, " ...")
    g_grid, L = line_profile_kernel(LAMP_PARAMS, :lamppost, :thin)

    println("Convolving (matrix path) ...")
    t_matrix = @elapsed F_matrix =
        convolve_reflection_matrix(R, log_lo, log_hi, g_grid, L; n_sub = 4)
    println("Convolving (FFT path) ...")
    convolve_reflection_fft(R, log_lo, log_hi, g_grid, L)  # warm-up (plan + compile)
    t_fft = @elapsed F_fft = convolve_reflection_fft(R, log_lo, log_hi, g_grid, L)

    rel_l2 = norm(F_fft - F_matrix) / max(norm(F_matrix), eps())
    println("matrix: ", round(t_matrix * 1e3, digits = 1), " ms")
    println("fft:    ", round(t_fft * 1e3, digits = 1), " ms (warm)")
    println("FFT vs matrix relative L2: ", rel_l2)
    assert_finite_nonneg("matrix flux", F_matrix)
    assert_finite_nonneg("fft flux", F_fft)
    rel_l2 <= 0.01 || error("FFT vs matrix disagree: relative L2 = $rel_l2 > 0.01")

    E = energy_midpoints(log_lo, log_hi)
    Fm = flux_per_keV(F_matrix, log_lo, log_hi)
    Ff = flux_per_keV(F_fft, log_lo, log_hi)
    denom = max.(Fm, 1e-12 * maximum(Fm))

    p_top = plot(
        E,
        positive_floor(Fm);
        xscale = :log10,
        yscale = :log10,
        ylabel = "Flux density (per keV)",
        label = "matrix ($(round(t_matrix * 1e3, digits = 1)) ms)",
        linewidth = 2,
        title = "Blur via matrix vs FFT (relative L2 = $(round(rel_l2, sigdigits = 2)))",
    )
    plot!(p_top, E, positive_floor(Ff);
        label = "fft ($(round(t_fft * 1e3, digits = 1)) ms)",
        linewidth = 2, linestyle = :dash)

    p_bottom = plot(
        E,
        (Ff .- Fm) ./ denom;
        xscale = :log10,
        xlabel = "Energy (keV)",
        ylabel = "(fft − matrix) / matrix",
        label = "",
        linewidth = 1,
    )
    hline!(p_bottom, [0.0]; color = :gray, linestyle = :dash, label = "")

    p = plot(p_top, p_bottom; layout = (2, 1), link = :x, size = (700, 600))
    save_plot(p, "ex4_convolution_methods.png")

    println("Done: 1 figure in ", OUT_DIR)
end

main()
