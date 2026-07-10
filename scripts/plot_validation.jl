#!/usr/bin/env julia
# Generate visual validation plots for the Julia-side pipeline.
#
# Usage (from repository root):
#   julia --project=. -e 'using Pkg; Pkg.add("Plots")'   # one-time
#   julia --project=. scripts/plot_validation.jl [path/to/xillverD-5.fits]

try
    using Plots
catch
    error(
        "Plots.jl is required. Run once:\n" *
        "  julia --project=. -e 'using Pkg; Pkg.add(\"Plots\")'",
    )
end

using GradusXSPEC

const DEFAULT_TABLE = "xillverD-5.fits"
const OUT_DIR = joinpath(@__DIR__, "..", "docs", "figures")

const GRADUS_PARAMS = (spin = 0.998, Eddington = 0.1, inc = 30.0, h = 3.0)
const REFL_PARAMS = (Gamma = 2.0, A_Fe = 1.0, logXi = 2.0, Dens = 17.0, Incl = 45.0)

function _integrate_g(g, f)
    area = 0.0
    for i in 1:(length(g) - 1)
        area += 0.5 * (f[i] + f[i + 1]) * (g[i + 1] - g[i])
    end
    return area
end

function narrow_kernel(g_grid; g0::Float64 = 1.0, sigma::Float64 = 0.001)
    g = collect(Float64.(g_grid))
    L = exp.(-0.5 .* ((g .- g0) ./ sigma).^2)
    return g, L ./ _integrate_g(g, L)
end

function energy_midpoints(energy_lo, energy_hi)
    return (Float64.(energy_lo) .+ Float64.(energy_hi)) ./ 2
end

function bin_widths(energy_lo, energy_hi)
    return Float64.(energy_hi) .- Float64.(energy_lo)
end

"""
    flux_per_keV(flux, energy_lo, energy_hi)

Convert XSPEC-style integrated bin fluxes to flux density for plotting only.
"""
function flux_per_keV(flux, energy_lo, energy_hi)
    widths = bin_widths(energy_lo, energy_hi)
    return Float64.(flux) ./ widths
end

function mask_energy_range(energy_lo, energy_hi, e_min::Float64, e_max::Float64)
    mid = energy_midpoints(energy_lo, energy_hi)
    return findall(e_min .<= mid .<= e_max)
end

function save_plot(plot, filename::String)
    path = joinpath(OUT_DIR, filename)
    mkpath(OUT_DIR)
    savefig(plot, path)
    println("Wrote ", path)
end

function plot_line_profile_kernel(gradus_params)
    g, L = line_profile_kernel(gradus_params)
    area = _integrate_g(g, L)
    p = plot(
        g,
        L;
        xlabel = "g = ν_obs / ν_em",
        ylabel = "L(g)",
        title = "Line-profile kernel (unit area = $(round(area, digits=4)))",
        legend = false,
        linewidth = 2,
    )
    vline!([1.0]; linestyle = :dash, color = :gray, label = "g = 1")
    save_plot(p, "01_line_profile_kernel.png")
end

function plot_reflection_spectrum(em_lo, em_hi, R)
    E = energy_midpoints(em_lo, em_hi)
    R_density = flux_per_keV(R, em_lo, em_hi)

    p_full = plot(
        E,
        R_density;
        xscale = :log10,
        yscale = :log10,
        xlabel = "Energy (keV)",
        ylabel = "Flux density (per keV)",
        title = "xillver reflection spectrum (full grid)",
        legend = false,
        linewidth = 1,
    )
    save_plot(p_full, "02_reflection_spectrum.png")

    idx = mask_energy_range(em_lo, em_hi, 3.0, 9.0)
    p_fe = plot(
        E[idx],
        R_density[idx];
        xlabel = "Energy (keV)",
        ylabel = "Flux density (per keV)",
        title = "xillver reflection (Fe Kα region, 3–9 keV)",
        legend = false,
        linewidth = 2,
    )
    save_plot(p_fe, "03_reflection_fe_region.png")
end

function plot_convolution_comparison(em_lo, em_hi, R, g_grid, L, title, filename)
    println("Convolving for plot: ", filename, " ...")
    F = convolve_reflection(R, em_lo, em_hi, g_grid, L; n_sub = 4)
    E = energy_midpoints(em_lo, em_hi)
    idx = mask_energy_range(em_lo, em_hi, 3.0, 9.0)
    R_plot = flux_per_keV(R, em_lo, em_hi)
    F_plot = flux_per_keV(F, em_lo, em_hi)

    p = plot(
        E[idx],
        R_plot[idx];
        label = "R(E) rest frame",
        xlabel = "Energy (keV)",
        ylabel = "Flux density (per keV)",
        title = title,
        linewidth = 2,
    )
    plot!(p, E[idx], F_plot[idx]; label = "F(E) convolved", linewidth = 2)
    save_plot(p, filename)
end

function plot_rebin_example(em_lo, em_hi, flux)
    n_dst = max(32, length(flux) ÷ 16)
    dst_edges = exp.(range(log(first(em_lo)), log(last(em_hi)), length = n_dst + 1))
    dst_lo = dst_edges[1:end-1]
    dst_hi = dst_edges[2:end]
    rebinned = rebin_flux(flux, em_lo, em_hi, dst_lo, dst_hi)

    E_src = energy_midpoints(em_lo, em_hi)
    E_dst = energy_midpoints(dst_lo, dst_hi)
    idx_src = mask_energy_range(em_lo, em_hi, 3.0, 9.0)
    idx_dst = mask_energy_range(dst_lo, dst_hi, 3.0, 9.0)
    flux_plot = flux_per_keV(flux, em_lo, em_hi)
    rebinned_plot = flux_per_keV(rebinned, dst_lo, dst_hi)

    p = plot(
        E_src[idx_src],
        flux_plot[idx_src];
        label = "table grid",
        xlabel = "Energy (keV)",
        ylabel = "Flux density (per keV)",
        title = "Flux-conserving rebin (3–9 keV)",
        linewidth = 1,
        xlims = (3.0, 9.0),
    )
    plot!(p, E_dst[idx_dst], rebinned_plot[idx_dst]; label = "rebinned ($(length(rebinned)) bins)", linewidth = 2)
    save_plot(p, "06_rebin_example.png")
end

function main()
    gradus_params = (
        Float64(GRADUS_PARAMS.spin),
        Float64(GRADUS_PARAMS.Eddington),
        Float64(GRADUS_PARAMS.inc),
        Float64(GRADUS_PARAMS.h),
    )
    refl_params = (
        Float64(REFL_PARAMS.Gamma),
        Float64(REFL_PARAMS.A_Fe),
        Float64(REFL_PARAMS.logXi),
        Float64(REFL_PARAMS.Dens),
        Float64(REFL_PARAMS.Incl),
    )

    table_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_TABLE
    isfile(table_path) || error("table model not found: $table_path")

    println("Gradus parameters: ", GRADUS_PARAMS)
    println("Reflection parameters: ", REFL_PARAMS)
    println("Output directory: ", OUT_DIR)

    println("Evaluating line-profile kernel (Gradus; may take ~1–2 min) ...")
    plot_line_profile_kernel(gradus_params)

    println("Loading xillver table ...")
    table = load_xspec_table(table_path)
    R = interpolate_table_spectrum(table, refl_params)
    em_lo, em_hi = table.energy_lo, table.energy_hi
    plot_reflection_spectrum(em_lo, em_hi, R)

    g_narrow, L_narrow = narrow_kernel(default_g_grid(; n = 512, g_min = 0.5, g_max = 1.5))
    plot_convolution_comparison(
        em_lo,
        em_hi,
        R,
        g_narrow,
        L_narrow,
        "Convolution with narrow kernel (3–9 keV)",
        "04_convolution_narrow_kernel.png",
    )

    g_gradus, L_gradus = line_profile_kernel(gradus_params)
    plot_convolution_comparison(
        em_lo,
        em_hi,
        R,
        g_gradus,
        L_gradus,
        "Convolution with Gradus lamppost kernel (3–9 keV)",
        "05_convolution_gradus_kernel.png",
    )

    println("Building rebin example from Gradus-blurred spectrum ...")
    F_gradus = convolve_reflection(R, em_lo, em_hi, g_gradus, L_gradus; n_sub = 4)
    plot_rebin_example(em_lo, em_hi, F_gradus)

    println("Done. See docs/visual_validation.md")
end

main()
