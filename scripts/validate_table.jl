#!/usr/bin/env julia
# Validate xillver table loading and interpolation.
#
# Usage (from repository root):
#   julia --project=. scripts/validate_table.jl [path/to/xillverD-5.fits]

using FITSIO

include(joinpath(@__DIR__, "..", "src", "table_model.jl"))

const DEFAULT_TABLE = "xillverD-5.fits"
const TEST_PARAMS = (
    Gamma = 2.0,
    A_Fe = 1.0,
    logXi = 2.0,
    Dens = 17.0,
    Incl = 45.0,
)

function main()
    table_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_TABLE
    println("Loading table: ", table_path)
    table = load_xspec_table(table_path)

    println("Model: ", table.name)
    println("Parameters: ", join(table.param_names, ", "))
    for (name, values) in zip(table.param_names, table.param_values)
        println("  $name: n=$(length(values)), range=$(extrema(values))")
    end
    println("Energy bins: ", length(table.energy_lo))
    println("Energy range (keV): ", extrema(table.energy_lo))

    params = (
        Float64(TEST_PARAMS.Gamma),
        Float64(TEST_PARAMS.A_Fe),
        Float64(TEST_PARAMS.logXi),
        Float64(TEST_PARAMS.Dens),
        Float64(TEST_PARAMS.Incl),
    )
    println("Interpolating at ", TEST_PARAMS)
    spectrum = interpolate_table_spectrum(table, params)
    println("Spectrum length: ", length(spectrum))
    println("Spectrum integral (sum): ", sum(spectrum))
    println("Spectrum max: ", maximum(spectrum))
    println("Done.")
end

main()
