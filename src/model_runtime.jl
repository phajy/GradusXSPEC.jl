struct ModelRuntime
    definition::XspecModelDefinition
    gradus_grids::Vector{Vector{Float64}}
    reflection_grids::Vector{Vector{Float64}}
end

function ModelRuntime(def::XspecModelDefinition)
    gradus_grids = [build_parameter_grid(spec) for spec in def.gradus_parameters]
    reflection_grids = [build_parameter_grid(spec) for spec in def.reflection_parameters]
    return ModelRuntime(def, gradus_grids, reflection_grids)
end

function n_gradus_params(rt::ModelRuntime)
    return length(rt.definition.gradus_parameters)
end

function n_reflection_params(rt::ModelRuntime)
    return length(rt.definition.reflection_parameters)
end

function n_physics_params(rt::ModelRuntime)
    return n_physics_params(rt.definition)
end

function physics_parameters(rt::ModelRuntime)
    return physics_parameters(rt.definition)
end

function _split_physics_params(rt::ModelRuntime, params::Tuple{Vararg{Float64, N}}) where {N}
    ng = n_gradus_params(rt)
    gradus = ntuple(i -> params[i], ng)
    reflection = ntuple(i -> params[ng + i], n_reflection_params(rt))
    return gradus, reflection
end

function _gradus_grid_point_params(rt::ModelRuntime, idx::Tuple{Vararg{Int, N}}) where {N}
    return ntuple(i -> rt.gradus_grids[i][idx[i]], N)
end

function _reflection_grid_point_params(rt::ModelRuntime, idx::Tuple{Vararg{Int, N}}) where {N}
    return ntuple(i -> rt.reflection_grids[i][idx[i]], N)
end

function _format_gradus_grid_point(rt::ModelRuntime, params::Tuple{Vararg{Float64, N}}) where {N}
    parts = String[]
    for (i, spec) in pairs(rt.definition.gradus_parameters)
        push!(parts, "$(spec.name)=$(params[i])")
    end
    return join(parts, ", ")
end

function _format_reflection_grid_point(rt::ModelRuntime, params::Tuple{Vararg{Float64, N}}) where {N}
    parts = String[]
    for (i, spec) in pairs(rt.definition.reflection_parameters)
        push!(parts, "$(spec.name)=$(params[i])")
    end
    return join(parts, ", ")
end

function _convolution_cache_signature(g_grid::AbstractVector{<:Real}, n_sub::Int)
    return UInt64(hash(g_grid, hash(n_sub)))
end

const MATRIX_CACHE_LOCK = ReentrantLock()
const CONVOLUTION_MATRIX_CACHE =
    Dict{Tuple{String, UInt64, Int, String, Tuple{Vararg{Int}}}, Matrix{Float32}}()
const MATRIX_CACHE_ORDER = Any[]
const MATRIX_CACHE_SIZES = Dict{Any, UInt64}()
const MATRIX_CACHE_HITS = Ref(0)
const MATRIX_CACHE_MISSES = Ref(0)

function _evict_one_matrix!()
    return _evict_one_from!(
        MATRIX_CACHE_LOCK,
        CONVOLUTION_MATRIX_CACHE,
        MATRIX_CACHE_ORDER,
        MATRIX_CACHE_SIZES,
    )
end

"""Apply a cached Float32 blur matrix to a Float64 reflection spectrum."""
function _apply_blur_matrix(M::Matrix{Float32}, R::AbstractVector{<:Real})
    return Float64.(M * Float32.(R))
end

function _get_or_compute_convolution_matrix(
    rt::ModelRuntime,
    table::XspecTableModel,
    table_path::AbstractString,
    gradus_idx::Tuple{Vararg{Int, N}};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
) where {N}
    sig = _convolution_cache_signature(g_grid, n_sub)
    key = (rt.definition.name, sig, n_sub, table_path, gradus_idx)

    cached = _bounded_cache_lookup!(
        MATRIX_CACHE_LOCK,
        CONVOLUTION_MATRIX_CACHE,
        MATRIX_CACHE_ORDER,
        key,
    )
    if cached !== nothing
        MATRIX_CACHE_HITS[] += 1
        return cached
    end

    MATRIX_CACHE_MISSES[] += 1
    L = _get_or_compute_line_kernel(rt, gradus_idx; g_grid = g_grid)
    M64 = build_convolution_matrix(
        table.energy_lo,
        table.energy_hi,
        table.energy_lo,
        table.energy_hi,
        g_grid,
        L;
        n_sub = n_sub,
    )
    M = Float32.(M64)
    return _bounded_cache_put!(
        MATRIX_CACHE_LOCK,
        CONVOLUTION_MATRIX_CACHE,
        MATRIX_CACHE_ORDER,
        MATRIX_CACHE_SIZES,
        key,
        M,
        _array_nbytes(M),
    )
end

function _interpolate_reflection_spectrum(
    rt::ModelRuntime,
    table::XspecTableModel,
    refl_corners::Dict{Tuple{Vararg{Int, N}}, Float64},
) where {N}
    n_bins = length(table.energy_lo)
    R = zeros(Float64, n_bins)
    for (refl_idx, weight) in refl_corners
        refl_params = _reflection_grid_point_params(rt, refl_idx)
        R_corner = interpolate_table_spectrum(table, refl_params)
        @inbounds for j in 1:n_bins
            R[j] += weight * R_corner[j]
        end
    end
    return R
end

function evaluate_spectrum_interpolated(
    rt::ModelRuntime,
    energy_edges::AbstractVector{<:Real},
    params::Tuple{Vararg{Float64, N}};
    table_path::AbstractString = DEFAULT_TABLE_PATH,
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
    verbose::Bool = false,
) where {N}
    table = get_table(table_path)
    resolved_path = _resolve_table_path(table_path)
    gradus_params, refl_params = _split_physics_params(rt, params)
    refl_corners = _multilinear_corners(
        Tuple(rt.reflection_grids),
        refl_params,
    )
    gradus_corners = _multilinear_corners(
        Tuple(rt.gradus_grids),
        gradus_params,
    )

    R_interp = _interpolate_reflection_spectrum(rt, table, refl_corners)

    used_identity = false
    if rt.definition.corona_variant == :gauss
        # The Gaussian diagnostic kernel is analytic and cheap, so evaluate it
        # at the exact requested Sigma instead of interpolating between grid
        # corners. Mixing an identity (narrow-σ) corner with a blurred corner
        # is not equivalent to blurring at the requested σ.
        #
        # Narrow Gaussians are a true identity: do not build a discrete δ(g-1)
        # kernel matrix (that path does not preserve counts and introduces an
        # ~E-dependent scale error versus atable{xillver}).
        if gaussian_is_identity(gradus_params[1]; g_grid = g_grid)
            used_identity = true
            convolved = copy(R_interp)
        else
            _, L = line_profile_kernel(
                gradus_params,
                rt.definition.corona_variant,
                rt.definition.disc_variant;
                g_grid = g_grid,
            )
            M = build_convolution_matrix(
                table.energy_lo,
                table.energy_hi,
                table.energy_lo,
                table.energy_hi,
                g_grid,
                L;
                n_sub = n_sub,
            )
            convolved = M * R_interp
        end
    else
        convolved = zeros(Float64, length(R_interp))
        for (gradus_idx, weight) in gradus_corners
            M = _get_or_compute_convolution_matrix(
                rt,
                table,
                resolved_path,
                gradus_idx;
                g_grid = g_grid,
                n_sub = n_sub,
            )
            blurred = _apply_blur_matrix(M, R_interp)
            @inbounds for i in eachindex(convolved)
                convolved[i] += weight * blurred[i]
            end
        end
    end

    output = _rebin_to_energy_edges(convolved, table.energy_lo, table.energy_hi, energy_edges)

    if verbose
        mode = used_identity ? "identity (no blur)" : "blur"
        println(
            "GradusXSPEC: $(rt.definition.name) $mode at ($(_format_gradus_grid_point(rt, gradus_params)); ",
            "$(_format_reflection_grid_point(rt, refl_params))) from ",
            "$(length(gradus_corners)) Gradus and $(length(refl_corners)) reflection corner(s); ",
            "matrix cache hits=$(MATRIX_CACHE_HITS[]), misses=$(MATRIX_CACHE_MISSES[])",
        )
    end

    _monitor_after_evaluation(rt, params, gradus_corners, refl_corners)

    return output
end

# Direct (non-interpolated) evaluation at the exact requested parameters.
# Used for validation against the grid-interpolated path; bypasses the
# parameter grids and convolution-matrix cache entirely.
function evaluate_spectrum(
    rt::ModelRuntime,
    energy_edges::AbstractVector{<:Real},
    params::Tuple{Vararg{Float64, N}};
    table_path::AbstractString = DEFAULT_TABLE_PATH,
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
    verbose::Bool = false,
) where {N}
    table = get_table(table_path)
    gradus_params, refl_params = _split_physics_params(rt, params)

    R = interpolate_table_spectrum(table, refl_params)

    if rt.definition.corona_variant == :gauss &&
       gaussian_is_identity(gradus_params[1]; g_grid = g_grid)
        convolved = copy(R)
    else
        _, L = line_profile_kernel(
            gradus_params,
            rt.definition.corona_variant,
            rt.definition.disc_variant;
            g_grid = g_grid,
        )
        convolved = convolve_reflection(
            R,
            table.energy_lo,
            table.energy_hi,
            g_grid,
            L;
            n_sub = n_sub,
        )
    end

    if verbose
        println(
            "GradusXSPEC: $(rt.definition.name) direct evaluation at ",
            "($(_format_gradus_grid_point(rt, gradus_params)); ",
            "$(_format_reflection_grid_point(rt, refl_params)))",
        )
    end

    return _rebin_to_energy_edges(convolved, table.energy_lo, table.energy_hi, energy_edges)
end

function build_model_runtimes(models = ALL_MODELS)
    return Dict(def.name => ModelRuntime(def) for def in models)
end

const MODEL_RUNTIMES = build_model_runtimes()

function get_model_runtime(name::AbstractString)
    rt = get(MODEL_RUNTIMES, name, nothing)
    rt === nothing && throw(ArgumentError("unknown model: $name"))
    return rt
end

# Backwards-compatible aliases for the original S&S model.
const GRADUS_PARAM_GRIDS = MODEL_RUNTIMES[LAMP_SS_MODEL.name].gradus_grids
const REFLECTION_PARAM_GRIDS = MODEL_RUNTIMES[LAMP_SS_MODEL.name].reflection_grids

function _gradus_grid_point_params(idx::NTuple{N_GRADUS_PARAMS, Int})
    return _gradus_grid_point_params(MODEL_RUNTIMES[LAMP_SS_MODEL.name], idx)
end

function _format_gradus_grid_point(params::NTuple{N_GRADUS_PARAMS, Float64})
    return _format_gradus_grid_point(MODEL_RUNTIMES[LAMP_SS_MODEL.name], params)
end

function _format_reflection_grid_point(params::NTuple{N_REFLECTION_PARAMS, Float64})
    return _format_reflection_grid_point(MODEL_RUNTIMES[LAMP_SS_MODEL.name], params)
end

function _split_physics_params(params::NTuple{N_XSPEC_FUNC_PARAMS, Float64})
    return _split_physics_params(MODEL_RUNTIMES[LAMP_SS_MODEL.name], params)
end

function evaluate_spectrum_interpolated(
    energy_edges::AbstractVector{<:Real},
    params::NTuple{N_XSPEC_FUNC_PARAMS, Float64};
    kwargs...,
)
    return evaluate_spectrum_interpolated(
        MODEL_RUNTIMES[LAMP_SS_MODEL.name],
        energy_edges,
        params;
        kwargs...,
    )
end

function evaluate_spectrum(
    energy_edges::AbstractVector{<:Real},
    params::NTuple{N_XSPEC_FUNC_PARAMS, Float64};
    kwargs...,
)
    return evaluate_spectrum(MODEL_RUNTIMES[LAMP_SS_MODEL.name], energy_edges, params; kwargs...)
end

function evaluate_spectrum(
    energy_edges::AbstractVector{<:Real},
    gradus_params::NTuple{4, Float64},
    refl_params::NTuple{5, Float64};
    kwargs...,
)
    params = (gradus_params..., refl_params...)
    return evaluate_spectrum(MODEL_RUNTIMES[LAMP_SS_MODEL.name], energy_edges, params; kwargs...)
end

function convolution_matrix_cache_size()
    return length(CONVOLUTION_MATRIX_CACHE)
end

function convolution_matrix_cache_stats()
    return MATRIX_CACHE_HITS[], MATRIX_CACHE_MISSES[]
end
