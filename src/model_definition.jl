const MODEL_NAME = "gradus"
const MODEL_TYPE = "add"
const MODEL_ENERGY_RANGE = (0.0, 1.0e20)
const MODEL_LANGUAGE_SYMBOL = "c_gradusjulia"

const MODEL_PARAMETERS = (
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

function model_dat_text()
    n_params = length(MODEL_PARAMETERS)
    lines = String[]
    push!(
        lines,
        "$(MODEL_NAME) $(n_params) $(MODEL_ENERGY_RANGE[1]) $(MODEL_ENERGY_RANGE[2]) $(MODEL_LANGUAGE_SYMBOL) $(MODEL_TYPE) 0 0",
    )
    for spec in MODEL_PARAMETERS
        unit = isempty(spec.unit) ? "\" \"" : spec.unit
        push!(
            lines,
            "$(spec.name) $(unit) $(spec.initial) $(spec.soft_min) $(spec.hard_min) $(spec.soft_max) $(spec.hard_max) $(spec.delta)",
        )
    end
    return join(lines, "\n") * "\n"
end
