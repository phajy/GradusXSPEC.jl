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

const PHYSICS_PARAM_GRIDS =
    ntuple(i -> build_parameter_grid(PHYSICS_PARAMETERS[i]), N_PHYSICS_PARAMS)
const REFLECTION_CACHE_LOCK = ReentrantLock()
const REFLECTION_SPECTRUM_CACHE =
    Dict{Tuple{UInt64, String, NTuple{N_PHYSICS_PARAMS, Int}}, Vector{Float64}}()
const REFLECTION_CACHE_HITS = Ref(0)
const REFLECTION_CACHE_MISSES = Ref(0)

function _physics_grid_bounds_and_weight(value::Float64, grid::Vector{Float64})
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

function _physics_interpolation_corners(params::NTuple{N_PHYSICS_PARAMS, Float64})
    bounds = ntuple(i -> _physics_grid_bounds_and_weight(params[i], PHYSICS_PARAM_GRIDS[i]), N_PHYSICS_PARAMS)
    corners = Dict{NTuple{N_PHYSICS_PARAMS, Int}, Float64}()
    for mask in 0:(UInt(1) << N_PHYSICS_PARAMS) - 1
        idx = ntuple(
            i -> ((mask >> (i - 1)) & UInt(1)) == UInt(1) ? bounds[i][2] : bounds[i][1],
            N_PHYSICS_PARAMS,
        )
        weight = 1.0
        for i in 1:N_PHYSICS_PARAMS
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

function _physics_grid_point_params(idx::NTuple{N_PHYSICS_PARAMS, Int})
    return ntuple(i -> PHYSICS_PARAM_GRIDS[i][idx[i]], N_PHYSICS_PARAMS)
end

function _format_physics_grid_point(params::NTuple{N_PHYSICS_PARAMS, Float64})
    parts = String[]
    for (i, spec) in pairs(PHYSICS_PARAMETERS)
        push!(parts, "$(spec.name)=$(params[i])")
    end
    return join(parts, ", ")
end

function _get_or_compute_convolved_spectrum(
    table::XspecTableModel,
    table_path::AbstractString,
    idx::NTuple{N_PHYSICS_PARAMS, Int};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
    n_sub::Int = 4,
)
    g_sig = hash(g_grid, hash(n_sub))
    key = (UInt64(g_sig), table_path, idx)
    grid_params = _physics_grid_point_params(idx)
    gradus_params = ntuple(i -> grid_params[i], N_GRADUS_PARAMS)
    refl_params = ntuple(i -> grid_params[N_GRADUS_PARAMS + i], N_REFLECTION_PARAMS)

    cached = lock(REFLECTION_CACHE_LOCK) do
        get(REFLECTION_SPECTRUM_CACHE, key, nothing)
    end
    if cached !== nothing
        REFLECTION_CACHE_HITS[] += 1
        return cached
    end

    REFLECTION_CACHE_MISSES[] += 1
    spec = _convolved_table_spectrum(
        table,
        gradus_params,
        refl_params;
        g_grid = g_grid,
        n_sub = n_sub,
    )
    lock(REFLECTION_CACHE_LOCK) do
        return get!(REFLECTION_SPECTRUM_CACHE, key, spec)
    end
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
    corners = _physics_interpolation_corners(params)
    n_bins = length(energy_edges) - 1
    output = zeros(Float64, n_bins)

    for (idx, weight) in corners
        convolved = _get_or_compute_convolved_spectrum(
            table,
            _resolve_table_path(table_path),
            idx;
            g_grid = g_grid,
            n_sub = n_sub,
        )
        rebinned = _rebin_to_energy_edges(convolved, table.energy_lo, table.energy_hi, energy_edges)
        @inbounds for i in eachindex(output)
            output[i] += weight * rebinned[i]
        end
    end

    if verbose
        println(
            "GradusXSPEC: reflection interpolation at ($(_format_physics_grid_point(params))) ",
            "from $(length(corners)) corner(s); cache hits=$(REFLECTION_CACHE_HITS[]), ",
            "misses=$(REFLECTION_CACHE_MISSES[])",
        )
    end

    return output
end
