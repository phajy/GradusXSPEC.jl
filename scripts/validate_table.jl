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

    println("Cross-checking SPECTRA rows against PARAMVAL ...")
    f = FITS(table_path)
    try
        pval = read(f[4], "PARAMVAL")
        intpspec = Float64.(read(f[4], "INTPSPEC"))
        nrows = size(pval, 2)
        mismatches = 0
        for r in 1:nrows
            idx = ntuple(i -> findfirst(==(Float64(pval[i, r])), table.param_values[i]), 5)
            from_table = vec(view(table.spectra, idx[5], idx[4], idx[3], idx[2], idx[1], :))
            direct = intpspec[:, r]
            if maximum(abs.(from_table .- direct)) > 1e-6 * max(maximum(abs.(direct)), 1.0)
                mismatches += 1
            end
        end
        println("PARAMVAL row mismatches: ", mismatches, " / ", nrows)
        mismatches == 0 || error("table indexing does not match PARAMVAL rows")
    finally
        close(f)
    end

    println("Done.")
end

main()
