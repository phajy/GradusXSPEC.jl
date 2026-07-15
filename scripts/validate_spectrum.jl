#!/usr/bin/env julia
# End-to-end Julia validation for the blurred reflection spectrum.
#
# Usage (from repository root):
#   julia --project=. scripts/validate_spectrum.jl [path/to/xillverD-5.fits]

using LinearAlgebra
using GradusXSPEC

const DEFAULT_TABLE = "xillverD-5.fits"

const GRADUS_PARAMS = (spin = 0.998, Eddington = 0.1, inc = 30.0, h = 3.0)
const REFL_PARAMS = (Gamma = 2.0, A_Fe = 1.0, logXi = 2.0, Dens = 17.0, Incl = 45.0)

function main()
    table_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_TABLE
    isfile(table_path) || error("table model not found: $table_path")

    gradus = (
        Float64(GRADUS_PARAMS.spin),
        Float64(GRADUS_PARAMS.Eddington),
        Float64(GRADUS_PARAMS.inc),
        Float64(GRADUS_PARAMS.h),
    )
    refl = (
        Float64(REFL_PARAMS.Gamma),
        Float64(REFL_PARAMS.A_Fe),
        Float64(REFL_PARAMS.logXi),
        Float64(REFL_PARAMS.Dens),
        Float64(REFL_PARAMS.Incl),
    )
    params = (gradus..., refl...)

    energy_edges = exp.(range(log(0.1), log(10.0), length = 101))
    println("Gradus parameters: ", GRADUS_PARAMS)
    println("Reflection parameters: ", REFL_PARAMS)
    println("Energy bins: ", length(energy_edges) - 1)

    println("Evaluating direct spectrum (Gradus + convolution; may take ~2 min) ...")
    direct = evaluate_spectrum(energy_edges, params; table_path = table_path)

    println("Evaluating interpolated spectrum (cached grid corners) ...")
    interpolated = evaluate_spectrum_interpolated(energy_edges, params; table_path = table_path)

    rel_diff = norm(direct - interpolated) / max(norm(interpolated), eps())
    println("Direct vs interpolated (relative L2): ", rel_diff)
    println("Direct flux sum: ", sum(direct))
    println("Interpolated flux sum: ", sum(interpolated))
    println("Direct flux max: ", maximum(direct))

    rel_diff < 0.05 || error("direct and interpolated spectra disagree")
    sum(direct) > 0 || error("spectrum flux is zero")
    println("Done.")
end

main()
