#!/usr/bin/env julia
# Validate line-profile kernel evaluation from Julia.
#
# Usage (from repository root):
#   julia --project=. scripts/validate_line_profile.jl

using GradusXSPEC

const TEST_PARAMS = (spin = 0.998, Eddington = 0.1, inc = 30.0, h = 3.0)

function _integrate_g(g, f)
    area = 0.0
    for i in 1:(length(g) - 1)
        area += 0.5 * (f[i] + f[i + 1]) * (g[i + 1] - g[i])
    end
    return area
end

function main()
    params = (
        Float64(TEST_PARAMS.spin),
        Float64(TEST_PARAMS.Eddington),
        Float64(TEST_PARAMS.inc),
        Float64(TEST_PARAMS.h),
    )
    println("Gradus parameters: ", TEST_PARAMS)

    g, kernel = line_profile_kernel(params)
    area = _integrate_g(g, kernel)
    println("Kernel grid: n=$(length(g)), g range=$(extrema(g))")
    println("Unit-area check (∫L dg): ", area)
    println("Kernel max: ", maximum(kernel))

    energies = exp.(range(log(0.1), log(10.0), length = 101))
    flux = line_profile_on_energy_edges(energies, params)
    println("Energy-edge profile: n=$(length(flux)) bins")
    println("Flux sum: ", sum(flux))
    println("Flux max: ", maximum(flux))
    println("Done.")
end

main()
