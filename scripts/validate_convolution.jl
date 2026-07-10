#!/usr/bin/env julia
# Validate reflection convolution and flux-conserving rebin.
#
# Usage (from repository root):
#   julia --project=. scripts/validate_convolution.jl [path/to/xillverD-5.fits]

using LinearAlgebra
using GradusXSPEC

const DEFAULT_TABLE = "xillverD-5.fits"
const TEST_REFL_PARAMS = (
    Gamma = 2.0,
    A_Fe = 1.0,
    logXi = 2.0,
    Dens = 17.0,
    Incl = 45.0,
)

function _integrate_g(g, f)
    area = 0.0
    for i in 1:(length(g) - 1)
        area += 0.5 * (f[i] + f[i + 1]) * (g[i + 1] - g[i])
    end
    return area
end

function narrow_kernel(g_grid; g0::Float64 = 1.0, sigma::Float64 = 0.01)
    g = collect(Float64.(g_grid))
    L = exp.(-0.5 .* ((g .- g0) ./ sigma).^2)
    area = _integrate_g(g, L)
    area > 0 || error("narrow kernel integral is zero")
    return L ./ area
end

function relative_l2(a, b)
    denom = max(norm(b), eps())
    return norm(a - b) / denom
end

function log_energy_edges(e_min::Float64, e_max::Float64, n_bins::Int)
    edges = exp.(range(log(e_min), log(e_max), length = n_bins + 1))
    return edges[1:end-1], edges[2:end]
end

function test_delta_kernel_on_log_grid()
    em_lo, em_hi = log_energy_edges(1.0, 20.0, 128)
    R = zeros(length(em_lo))
    R[64] = 1.0
    g_grid = default_g_grid(; n = 512, g_min = 0.5, g_max = 1.5)
    L = narrow_kernel(g_grid; g0 = 1.0, sigma = 0.001)
    F = convolve_reflection(R, em_lo, em_hi, g_grid, L; n_sub = 16)
    shape_err = relative_l2(F ./ max(sum(F), eps()), R ./ max(sum(R), eps()))
    println("Synthetic delta-kernel shape error (relative L2): ", shape_err)
    println("Synthetic peak bin shift: ", argmax(F) - argmax(R))
    return shape_err < 0.06 && argmax(F) == argmax(R)
end

function test_matrix_matches_direct(em_lo, em_hi, R, g_grid, L)
    M = build_convolution_matrix(em_lo, em_hi, em_lo, em_hi, g_grid, L; n_sub = 4)
    F_mat = M * R
    F_dir = convolve_reflection(R, em_lo, em_hi, g_grid, L; n_sub = 4)
    err = relative_l2(F_mat, F_dir)
    println("Matrix vs direct convolution error (relative L2): ", err)
    return err < 1e-14
end

function test_rebin_conserves_flux(em_lo, em_hi, flux)
    n_dst = max(1, length(flux) ÷ 4)
    dst_edges = exp.(range(log(first(em_lo)), log(last(em_hi)), length = n_dst + 1))
    dst_lo = dst_edges[1:end-1]
    dst_hi = dst_edges[2:end]
    rebinned = rebin_flux(flux, em_lo, em_hi, dst_lo, dst_hi)
    src_total = sum(flux)
    dst_total = sum(rebinned)
    rel_diff = abs(dst_total - src_total) / max(src_total, eps())
    println("Rebin flux conservation (relative difference): ", rel_diff)
    return rel_diff < 1e-12
end

function main()
    ok_synthetic = test_delta_kernel_on_log_grid()

    table_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_TABLE
    isfile(table_path) || error("table model not found: $table_path")

    println("Loading table: ", table_path)
    table = load_xspec_table(table_path)
    refl_params = (
        Float64(TEST_REFL_PARAMS.Gamma),
        Float64(TEST_REFL_PARAMS.A_Fe),
        Float64(TEST_REFL_PARAMS.logXi),
        Float64(TEST_REFL_PARAMS.Dens),
        Float64(TEST_REFL_PARAMS.Incl),
    )
    R = interpolate_table_spectrum(table, refl_params)
    em_lo, em_hi = table.energy_lo, table.energy_hi
    println("Reflection spectrum bins: ", length(R))

    g_grid = default_g_grid(; n = 128)
    L = narrow_kernel(g_grid; g0 = 1.0, sigma = 0.05)
    ok_matrix = test_matrix_matches_direct(em_lo, em_hi, R, g_grid, L)
    ok_rebin = test_rebin_conserves_flux(em_lo, em_hi, R)

    if ok_synthetic && ok_matrix && ok_rebin
        println("All convolution checks passed.")
    else
        error("Convolution validation failed.")
    end
end

main()
