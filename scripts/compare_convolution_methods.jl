#!/usr/bin/env julia
# Compare matrix vs FFT reflection convolution on the same inputs.
#
# Usage (from repository root):
#   julia --project=. scripts/compare_convolution_methods.jl [path/to/xillverD-5.fits]

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

function narrow_kernel(g_grid; g0::Float64 = 1.0, sigma::Float64 = 0.05)
    g = collect(Float64.(g_grid))
    L = exp.(-0.5 .* ((g .- g0) ./ sigma) .^ 2)
    area = _integrate_g(g, L)
    area > 0 || error("narrow kernel integral is zero")
    return L ./ area
end

function relative_l2(a, b)
    return norm(a - b) / max(norm(b), eps())
end

function compare_on(R, em_lo, em_hi, g_grid, L; label::AbstractString, n_sub::Int = 4, tol::Float64 = 0.01)
    println("--- ", label, " ---")
    Fm = convolve_reflection_matrix(R, em_lo, em_hi, g_grid, L; n_sub = n_sub)
    Ff = convolve_reflection_fft(R, em_lo, em_hi, g_grid, L)
    rel = relative_l2(Ff, Fm)
    println("  sum(R)=", sum(R), "  sum(matrix)=", sum(Fm), "  sum(fft)=", sum(Ff))
    println("  flux ratios: matrix/R=", sum(Fm) / sum(R), "  fft/R=", sum(Ff) / sum(R))
    println("  FFT vs matrix relative L2: ", rel)
    rel <= tol || error("FFT vs matrix disagree (relative L2=$rel > $tol) for $label")
    return rel
end

function main()
    println("Default convolution_method(): ", convolution_method())

    # Synthetic spike on a log grid (cleanest comparison).
    edges = exp.(range(log(0.5), log(20.0), length = 257))
    em_lo, em_hi = edges[1:end-1], edges[2:end]
    R = zeros(length(em_lo))
    R[128] = 1.0
    g_grid = default_g_grid(; n = 512, g_min = 0.5, g_max = 1.5)
    L = narrow_kernel(g_grid; sigma = 0.05)
    compare_on(R, em_lo, em_hi, g_grid, L; label = "synthetic spike, sigma=0.05")

    table_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_TABLE
    isfile(table_path) || error("table model not found: $table_path")
    println("Loading table: ", table_path)
    table = load_xspec_table(table_path)
    refl = (
        Float64(TEST_REFL_PARAMS.Gamma),
        Float64(TEST_REFL_PARAMS.A_Fe),
        Float64(TEST_REFL_PARAMS.logXi),
        Float64(TEST_REFL_PARAMS.Dens),
        Float64(TEST_REFL_PARAMS.Incl),
    )
    Rtab = interpolate_table_spectrum(table, refl)
    # Compare on a uniform log grid (FFT's native setting), not the hybrid blur grid.
    log_lo, log_hi = log_energy_bin_edges(0.5, 20.0, 512)
    R_log = rebin_flux(Rtab, table.energy_lo, table.energy_hi, log_lo, log_hi)
    L2 = narrow_kernel(g_grid; sigma = 0.05)
    compare_on(
        R_log,
        log_lo,
        log_hi,
        g_grid,
        L2;
        label = "xillver on uniform log grid 0.5–20 keV, sigma=0.05",
        tol = 0.01,
    )

    println("All FFT vs matrix comparisons passed.")
end

main()
