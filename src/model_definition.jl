const PACKAGE_NAME = "gradusxspec"
const MODEL_TYPE = "add"
const MODEL_ENERGY_RANGE = (0.0, 1.0e20)
const DEFAULT_TABLE_PATH = "xillverD-5.fits"
const DEFAULT_MONITOR_PATH = joinpath(dirname(@__DIR__), "gradusxspec_monitor.txt")

const LAMP_SS_GRADUS_PARAMETERS = (
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
        # Cap at 65°: Gradus transfer functions currently fail for some
        # (spin, h) at inc ≥ 70° with maxrₑ = 400 (see reproduce_transfer_offset.jl).
        soft_max = 65.0,
        hard_max = 65.0,
        delta = 5.0,
        interpolation = :linear,
    ),
    (
        name = "h",
        unit = "r_g",
        initial = 3.0,
        soft_min = 2.0,
        hard_min = 2.0,
        soft_max = 20.0,
        hard_max = 20.0,
        delta = 0.25,
        interpolation = :linear,
    ),
)

const LAMP_THIN_GRADUS_PARAMETERS = (
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
        name = "inc",
        unit = "degrees",
        initial = 30.0,
        soft_min = 5.0,
        hard_min = 5.0,
        # Cap at 65°: Gradus transfer functions currently fail for some
        # (spin, h) at inc ≥ 70° with maxrₑ = 400 (see reproduce_transfer_offset.jl).
        soft_max = 65.0,
        hard_max = 65.0,
        delta = 5.0,
        interpolation = :linear,
    ),
    (
        name = "h",
        unit = "r_g",
        initial = 3.0,
        soft_min = 2.0,
        hard_min = 2.0,
        soft_max = 20.0,
        hard_max = 20.0,
        delta = 0.25,
        interpolation = :linear,
    ),
)

# Temporary diagnostic model: Gaussian blur in g = E_obs/E_em.
# With Sigma ≪ 1 the kernel is nearly a delta at g=1, so the output is
# essentially the interpolated xillver table (no relativistic blurring).
const TEST_GAUSS_PARAMETERS = (
    (
        name = "Sigma",
        unit = "",
        initial = 0.0001,
        soft_min = 0.0001,
        hard_min = 0.0001,
        soft_max = 0.2,
        hard_max = 0.2,
        delta = 0.01,
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

struct XspecModelDefinition
    name::String
    language_symbol::String
    cc_name::String
    gradus_parameters::Tuple
    reflection_parameters::Tuple
    disc_variant::Symbol
end

function physics_parameters(def::XspecModelDefinition)
    return (def.gradus_parameters..., def.reflection_parameters...)
end

function n_gradus_params(def::XspecModelDefinition)
    return length(def.gradus_parameters)
end

function n_reflection_params(def::XspecModelDefinition)
    return length(def.reflection_parameters)
end

function n_physics_params(def::XspecModelDefinition)
    return length(physics_parameters(def))
end

function n_model_dat_params(def::XspecModelDefinition)
    return n_physics_params(def) + 1
end

const LAMP_SS_MODEL = XspecModelDefinition(
    "gradus_lamp_ss",
    "c_graduslampsjulia",
    "graduslampsjxspec",
    LAMP_SS_GRADUS_PARAMETERS,
    REFLECTION_PARAMETERS,
    :ss,
)

const LAMP_THIN_MODEL = XspecModelDefinition(
    "gradus_lamp_thin",
    "c_graduslampthinjulia",
    "graduslampthinxspec",
    LAMP_THIN_GRADUS_PARAMETERS,
    REFLECTION_PARAMETERS,
    :thin,
)

const TEST_GAUSS_MODEL = XspecModelDefinition(
    "test_gauss",
    "c_testgaussjulia",
    "testgaussxspec",
    TEST_GAUSS_PARAMETERS,
    REFLECTION_PARAMETERS,
    :gauss,
)

const ALL_MODELS = (LAMP_SS_MODEL, LAMP_THIN_MODEL, TEST_GAUSS_MODEL)

# Backwards-compatible aliases for the original lamppost + S&S model.
const MODEL_NAME = LAMP_SS_MODEL.name
const GRADUS_PARAMETERS = LAMP_SS_GRADUS_PARAMETERS
const PHYSICS_PARAMETERS = physics_parameters(LAMP_SS_MODEL)
const N_GRADUS_PARAMS = n_gradus_params(LAMP_SS_MODEL)
const N_REFLECTION_PARAMS = n_reflection_params(LAMP_SS_MODEL)
const N_PHYSICS_PARAMS = n_physics_params(LAMP_SS_MODEL)
const N_MODEL_DAT_PARAMS = n_model_dat_params(LAMP_SS_MODEL)
const N_XSPEC_FUNC_PARAMS = N_PHYSICS_PARAMS
const MODEL_PARAMETERS = PHYSICS_PARAMETERS

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

function _model_dat_block(def::XspecModelDefinition; table_path::AbstractString = DEFAULT_TABLE_PATH)
    params = physics_parameters(def)
    n_params = n_model_dat_params(def)
    lines = String[
        "$(def.name) $(n_params) $(MODEL_ENERGY_RANGE[1]) $(MODEL_ENERGY_RANGE[2]) $(def.language_symbol) $(MODEL_TYPE) 0 $(table_path)",
    ]
    for spec in params
        unit = isempty(spec.unit) ? "\" \"" : spec.unit
        push!(
            lines,
            "$(spec.name) $(unit) $(spec.initial) $(spec.soft_min) $(spec.hard_min) $(spec.soft_max) $(spec.hard_max) $(spec.delta)",
        )
    end
    push!(lines, "*redshift \" \" 0.0")
    return join(lines, "\n")
end

function model_dat_text(; table_path::AbstractString = DEFAULT_TABLE_PATH, models = ALL_MODELS)
    # XSPEC requires a blank line between model stanzas: updateComponentList
    # only treats a line as a new model header when the previous line is blank.
    # Without that separator, later models are never registered and lmod aborts
    # with XSModelFunction::NoSuchComponent when wiring function pointers.
    blocks = [_model_dat_block(def; table_path = table_path) for def in models]
    return join(blocks, "\n\n") * "\n"
end
