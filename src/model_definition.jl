const MODEL_NAME = "gradus"
const MODEL_TYPE = "add"
const MODEL_ENERGY_RANGE = (0.0, 1.0e20)
const MODEL_LANGUAGE_SYMBOL = "c_gradusjulia"
const DEFAULT_TABLE_PATH = "xillverD-5.fits"

const GRADUS_PARAMETERS = (
    (
        name = "spin",
        unit = "",
        initial = 0.998,
        soft_min = 0.0,
        hard_min = 0.0,
        soft_max = 0.998,
        hard_max = 0.998,
        delta = 0.05,
        interpolation = :linear,
    ),
    (
        name = "Eddington",
        unit = "",
        initial = 0.1,
        soft_min = 0.0,
        hard_min = 0.0,
        soft_max = 1.0,
        hard_max = 1.0,
        delta = 0.1,
        interpolation = :linear,
    ),
    (
        name = "inc",
        unit = "degrees",
        initial = 30.0,
        soft_min = 5.0,
        hard_min = 5.0,
        soft_max = 85.0,
        hard_max = 85.0,
        delta = 5.0,
        interpolation = :linear,
    ),
    (
        name = "h",
        unit = "r_g",
        initial = 3.0,
        soft_min = 1.0,
        hard_min = 1.0,
        soft_max = 20.0,
        hard_max = 20.0,
        delta = 0.25,
        interpolation = :linear,
    ),
)

const REFLECTION_PARAMETERS = (
    (
        name = "Refl_Gamma",
        unit = "",
        initial = 2.0,
        soft_min = 1.2,
        hard_min = 1.2,
        soft_max = 3.6,
        hard_max = 3.6,
        delta = 0.01,
        interpolation = :linear,
    ),
    (
        name = "Refl_A_Fe",
        unit = "",
        initial = 1.0,
        soft_min = 0.5,
        hard_min = 0.5,
        soft_max = 20.0,
        hard_max = 20.0,
        delta = 0.01,
        interpolation = :linear,
    ),
    (
        name = "Refl_logXi",
        unit = "",
        initial = 2.0,
        soft_min = 0.0,
        hard_min = 0.0,
        soft_max = 4.69897,
        hard_max = 4.69897,
        delta = 0.02,
        interpolation = :linear,
    ),
    (
        name = "Refl_Dens",
        unit = "",
        initial = 17.0,
        soft_min = 15.0,
        hard_min = 15.0,
        soft_max = 19.0,
        hard_max = 19.0,
        delta = 0.1,
        interpolation = :linear,
    ),
    (
        name = "Refl_Incl",
        unit = "degrees",
        initial = 45.0,
        soft_min = 18.1949,
        hard_min = 18.1949,
        soft_max = 87.134,
        hard_max = 87.134,
        delta = 0.45,
        interpolation = :linear,
    ),
)

# Combined parameter tuple used once reflection convolution is wired into XSPEC.
const MODEL_PARAMETERS = (GRADUS_PARAMETERS..., REFLECTION_PARAMETERS...)

const N_GRADUS_PARAMS = length(GRADUS_PARAMETERS)
const N_REFLECTION_PARAMS = length(REFLECTION_PARAMETERS)

function build_parameter_grid(spec)
    minv = Float64(spec.hard_min)
    maxv = Float64(spec.hard_max)
    stepv = Float64(spec.delta)
    if stepv <= 0
        return [minv, maxv]
    end
    values = collect(minv:stepv:maxv)
    if isempty(values)
        values = [minv]
    end
    if values[end] < maxv
        push!(values, maxv)
    elseif values[end] > maxv
        values[end] = maxv
    end
    return values
end

function model_dat_text(; table_path::AbstractString = DEFAULT_TABLE_PATH, include_reflection::Bool = false)
    params = include_reflection ? MODEL_PARAMETERS : GRADUS_PARAMETERS
    n_params = length(params)
    init_string = include_reflection ? table_path : "0"
    lines = String[]
    push!(
        lines,
        "$(MODEL_NAME) $(n_params) $(MODEL_ENERGY_RANGE[1]) $(MODEL_ENERGY_RANGE[2]) $(MODEL_LANGUAGE_SYMBOL) $(MODEL_TYPE) 0 $(init_string)",
    )
    for spec in params
        unit = isempty(spec.unit) ? "\" \"" : spec.unit
        push!(
            lines,
            "$(spec.name) $(unit) $(spec.initial) $(spec.soft_min) $(spec.hard_min) $(spec.soft_max) $(spec.hard_max) $(spec.delta)",
        )
    end
    return join(lines, "\n") * "\n"
end
