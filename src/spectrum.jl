"""
    InitConfig

Parsed XSPEC model init string. The first token is the reflection table path;
`verbose` may appear as an additional token.
"""
struct InitConfig
    table_path::String
    verbose::Bool
end

function parse_init_string(init::AbstractString)
    tokens = filter(!isempty, split(strip(init)))
    table_path = DEFAULT_TABLE_PATH
    verbose = false
    for token in tokens
        if lowercase(token) == "verbose"
            verbose = true
        elseif token != "0"
            table_path = token
        end
    end
    return InitConfig(table_path, verbose)
end

function parse_init_string(init::Ptr{Cchar})
    init == C_NULL && return InitConfig(DEFAULT_TABLE_PATH, false)
    return parse_init_string(strip(unsafe_string(init)))
end

const TABLE_LOCK = ReentrantLock()
const LOADED_TABLE = Ref{Union{XspecTableModel,Nothing}}(nothing)
const LOADED_TABLE_PATH = Ref("")

function get_table(table_path::AbstractString = DEFAULT_TABLE_PATH)
    resolved = _resolve_table_path(table_path)
    lock(TABLE_LOCK) do
        if LOADED_TABLE[] === nothing || LOADED_TABLE_PATH[] != resolved
            LOADED_TABLE[] = load_xspec_table(resolved)
            LOADED_TABLE_PATH[] = resolved
        end
        return LOADED_TABLE[]
    end
end

function _energy_signature(energies::AbstractVector{<:Real})
    h = UInt(0)
    @inbounds for e in energies
        h = hash(Float64(e), h)
    end
    return UInt64(h)
end

function _split_physics_params(params::NTuple{N_XSPEC_FUNC_PARAMS, Float64})
    gradus = ntuple(i -> params[i], N_GRADUS_PARAMS)
    reflection = ntuple(i -> params[N_GRADUS_PARAMS + i], N_REFLECTION_PARAMS)
    return gradus, reflection
end

function _convolved_table_spectrum(
    table::XspecTableModel,
    gradus_params::NTuple{4, Float64},
    refl_params::NTuple{5, Float64};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
)
    R = interpolate_table_spectrum(table, refl_params)
    _, L = line_profile_kernel(gradus_params; g_grid = g_grid)
    return convolve_reflection(
        R,
        table.energy_lo,
        table.energy_hi,
        g_grid,
        L;
        n_sub = n_sub,
    )
end

function _rebin_to_energy_edges(
    flux::AbstractVector{<:Real},
    src_lo::AbstractVector{<:Real},
    src_hi::AbstractVector{<:Real},
    energy_edges::AbstractVector{<:Real},
)
    length(energy_edges) >= 2 || throw(ArgumentError("energy_edges must contain at least two values"))
    dst_lo = Float64.(energy_edges[1:(end - 1)])
    dst_hi = Float64.(energy_edges[2:end])
    return rebin_flux(flux, src_lo, src_hi, dst_lo, dst_hi)
end

"""
    evaluate_spectrum(energy_edges, params; table_path=DEFAULT_TABLE_PATH, g_grid=default_g_grid())

Evaluate the blurred reflection spectrum for the nine physics parameters
`(spin, Eddington, inc, h, Refl_Gamma, Refl_A_Fe, Refl_logXi, Refl_Dens,
Refl_Incl)`. Returns bin-integrated flux on `energy_edges` in the rest frame;
XSPEC applies `norm` and `*redshift` externally for the local model.
"""
function evaluate_spectrum(
    energy_edges::AbstractVector{<:Real},
    params::NTuple{N_XSPEC_FUNC_PARAMS, Float64};
    table_path::AbstractString = DEFAULT_TABLE_PATH,
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
)
    gradus_params, refl_params = _split_physics_params(params)
    table = get_table(table_path)
    convolved = _convolved_table_spectrum(
        table,
        gradus_params,
        refl_params;
        g_grid = g_grid,
        n_sub = n_sub,
    )
    return _rebin_to_energy_edges(convolved, table.energy_lo, table.energy_hi, energy_edges)
end

"""
    evaluate_spectrum(energy_edges, gradus_params, refl_params; kwargs...)

Convenience overload with separate Gradus and reflection parameter tuples.
"""
function evaluate_spectrum(
    energy_edges::AbstractVector{<:Real},
    gradus_params::NTuple{4, Float64},
    refl_params::NTuple{5, Float64};
    kwargs...,
)
    params = (gradus_params..., refl_params...)
    return evaluate_spectrum(energy_edges, params; kwargs...)
end

const GRADUS_PARAM_GRIDS =
    ntuple(i -> build_parameter_grid(GRADUS_PARAMETERS[i]), N_GRADUS_PARAMS)
const REFLECTION_PARAM_GRIDS =
    ntuple(i -> build_parameter_grid(REFLECTION_PARAMETERS[i]), N_REFLECTION_PARAMS)

function _param_grid_bounds_and_weight(value::Float64, grid::Vector{Float64})
    clamped = clamp(value, first(grid), last(grid))
    hi = searchsortedfirst(grid, clamped)
    if hi <= 1
        return (1, 1, 0.0)
    elseif hi > length(grid)
        idx = length(grid)
        return (idx, idx, 0.0)
    elseif grid[hi] == clamped
        return (hi, hi, 0.0)
    else
        lo = hi - 1
        θ = (clamped - grid[lo]) / (grid[hi] - grid[lo])
        return (lo, hi, θ)
    end
end

function _multilinear_corners(
    param_grids::NTuple{N, Vector{Float64}},
    params::NTuple{N, Float64},
) where {N}
    bounds = ntuple(i -> _param_grid_bounds_and_weight(params[i], param_grids[i]), N)
    corners = Dict{NTuple{N, Int}, Float64}()
    for mask in 0:(UInt(1) << N) - 1
        idx = ntuple(
            i -> ((mask >> (i - 1)) & UInt(1)) == UInt(1) ? bounds[i][2] : bounds[i][1],
            N,
        )
        weight = 1.0
        for i in 1:N
            θ = bounds[i][3]
            if ((mask >> (i - 1)) & UInt(1)) == UInt(1)
                weight *= θ
            else
                weight *= (1 - θ)
            end
        end
        if weight > 0
            corners[idx] = get(corners, idx, 0.0) + weight
        end
    end
    return corners
end

function _gradus_grid_point_params(idx::NTuple{N_GRADUS_PARAMS, Int})
    return ntuple(i -> GRADUS_PARAM_GRIDS[i][idx[i]], N_GRADUS_PARAMS)
end

function _reflection_grid_point_params(idx::NTuple{N_REFLECTION_PARAMS, Int})
    return ntuple(i -> REFLECTION_PARAM_GRIDS[i][idx[i]], N_REFLECTION_PARAMS)
end

function _format_gradus_grid_point(params::NTuple{N_GRADUS_PARAMS, Float64})
    parts = String[]
    for (i, spec) in pairs(GRADUS_PARAMETERS)
        push!(parts, "$(spec.name)=$(params[i])")
    end
    return join(parts, ", ")
end

function _format_reflection_grid_point(params::NTuple{N_REFLECTION_PARAMS, Float64})
    parts = String[]
    for (i, spec) in pairs(REFLECTION_PARAMETERS)
        push!(parts, "$(spec.name)=$(params[i])")
    end
    return join(parts, ", ")
end

function _convolution_cache_signature(g_grid::AbstractVector{<:Real}, n_sub::Int)
    return UInt64(hash(g_grid, hash(n_sub)))
end

const MATRIX_CACHE_LOCK = ReentrantLock()
const CONVOLUTION_MATRIX_CACHE =
    Dict{Tuple{UInt64, Int, String, NTuple{N_GRADUS_PARAMS, Int}}, Matrix{Float64}}()
const MATRIX_CACHE_HITS = Ref(0)
const MATRIX_CACHE_MISSES = Ref(0)

function _get_or_compute_convolution_matrix(
    table::XspecTableModel,
    table_path::AbstractString,
    gradus_idx::NTuple{N_GRADUS_PARAMS, Int};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
)
    sig = _convolution_cache_signature(g_grid, n_sub)
    key = (sig, n_sub, table_path, gradus_idx)

    cached = lock(MATRIX_CACHE_LOCK) do
        get(CONVOLUTION_MATRIX_CACHE, key, nothing)
    end
    if cached !== nothing
        MATRIX_CACHE_HITS[] += 1
        return cached
    end

    MATRIX_CACHE_MISSES[] += 1
    gradus_params = _gradus_grid_point_params(gradus_idx)
    _, L = line_profile_kernel(gradus_params; g_grid = g_grid)
    M = build_convolution_matrix(
        table.energy_lo,
        table.energy_hi,
        table.energy_lo,
        table.energy_hi,
        g_grid,
        L;
        n_sub = n_sub,
    )
    lock(MATRIX_CACHE_LOCK) do
        return get!(CONVOLUTION_MATRIX_CACHE, key, M)
    end
end

function _interpolate_reflection_spectrum(
    table::XspecTableModel,
    refl_corners::Dict{NTuple{N_REFLECTION_PARAMS, Int}, Float64},
)
    n_bins = length(table.energy_lo)
    R = zeros(Float64, n_bins)
    for (refl_idx, weight) in refl_corners
        refl_params = _reflection_grid_point_params(refl_idx)
        R_corner = interpolate_table_spectrum(table, refl_params)
        @inbounds for j in 1:n_bins
            R[j] += weight * R_corner[j]
        end
    end
    return R
end

function evaluate_spectrum_interpolated(
    energy_edges::AbstractVector{<:Real},
    params::NTuple{N_XSPEC_FUNC_PARAMS, Float64};
    table_path::AbstractString = DEFAULT_TABLE_PATH,
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
    verbose::Bool = false,
)
    table = get_table(table_path)
    resolved_path = _resolve_table_path(table_path)
    gradus_params, refl_params = _split_physics_params(params)
    refl_corners = _multilinear_corners(REFLECTION_PARAM_GRIDS, refl_params)
    gradus_corners = _multilinear_corners(GRADUS_PARAM_GRIDS, gradus_params)

    R_interp = _interpolate_reflection_spectrum(table, refl_corners)
    convolved = zeros(Float64, length(R_interp))
    for (gradus_idx, weight) in gradus_corners
        M = _get_or_compute_convolution_matrix(
            table,
            resolved_path,
            gradus_idx;
            g_grid = g_grid,
            n_sub = n_sub,
        )
        blurred = M * R_interp
        @inbounds for i in eachindex(convolved)
            convolved[i] += weight * blurred[i]
        end
    end

    output = _rebin_to_energy_edges(convolved, table.energy_lo, table.energy_hi, energy_edges)

    if verbose
        println(
            "GradusXSPEC: reflection at ($(_format_gradus_grid_point(gradus_params)); ",
            "$(_format_reflection_grid_point(refl_params))) from ",
            "$(length(gradus_corners)) Gradus and $(length(refl_corners)) reflection corner(s); ",
            "matrix cache hits=$(MATRIX_CACHE_HITS[]), misses=$(MATRIX_CACHE_MISSES[])",
        )
    end

    return output
end
