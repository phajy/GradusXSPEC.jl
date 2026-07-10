module GradusXSPEC

using Base: @ccallable
using Gradus
include("model_definition.jl")
include("line_profile.jl")
include("table_model.jl")
include("convolution.jl")

export default_g_grid,
    g_midpoints_from_energy_edges,
    line_profile_kernel,
    line_profile_on_energy_edges,
    E_REST_KEV,
    load_xspec_table,
    interpolate_table_spectrum,
    interpolate_line_kernel,
    build_convolution_matrix,
    convolve_reflection,
    rebin_flux

const N_MODEL_PARAMS = length(GRADUS_PARAMETERS)
const PARAM_GRIDS = ntuple(i -> build_parameter_grid(GRADUS_PARAMETERS[i]), N_MODEL_PARAMS)
const CACHE_LOCK = ReentrantLock()
const SPECTRUM_CACHE = Dict{Tuple{UInt64,NTuple{N_MODEL_PARAMS,Int}},Vector{Float64}}()
const VERBOSE = Ref(false)
const CACHE_HITS = Ref(0)
const CACHE_MISSES = Ref(0)

function _verbose_enabled()
    return VERBOSE[] || lowercase(get(ENV, "GRADUSXSPEC_VERBOSE", "0")) in ("1", "true", "yes")
end

function _enable_verbose_from_init(init::Ptr{Cchar})
    init == C_NULL && return
    s = strip(unsafe_string(init))
    isempty(s) && return
    if occursin("verbose", lowercase(s))
        VERBOSE[] = true
    end
end

function _format_grid_point(params::NTuple{N_MODEL_PARAMS,Float64})
    parts = String[]
    for (i, spec) in pairs(GRADUS_PARAMETERS)
        push!(parts, "$(spec.name)=$(params[i])")
    end
    return join(parts, ", ")
end

function _energy_signature(energies::AbstractVector{<:Real})
    h = UInt(0)
    @inbounds for e in energies
        h = hash(Float64(e), h)
    end
    return UInt64(h)
end

# Core physics entry point: evaluates the Gradus line profile at one parameter-grid
# point. Grid caching and interpolation in `_interpolated_spectrum` call this.
function _evaluate_model_on_grid(params::NTuple{N_MODEL_PARAMS,Float64}, energies::AbstractVector{<:Real})
    return line_profile_on_energy_edges(energies, params)
end

function _grid_bounds_and_weight(value::Float64, grid::Vector{Float64})
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

function _interpolation_corners(params::NTuple{N_MODEL_PARAMS,Float64})
    bounds = ntuple(i -> _grid_bounds_and_weight(params[i], PARAM_GRIDS[i]), N_MODEL_PARAMS)
    corners = Dict{NTuple{N_MODEL_PARAMS,Int},Float64}()
    for mask in 0:(UInt(1) << N_MODEL_PARAMS) - 1
        idx = ntuple(i -> ((mask >> (i - 1)) & UInt(1)) == UInt(1) ? bounds[i][2] : bounds[i][1], N_MODEL_PARAMS)
        w = 1.0
        for i in 1:N_MODEL_PARAMS
            θ = bounds[i][3]
            if ((mask >> (i - 1)) & UInt(1)) == UInt(1)
                w *= θ
            else
                w *= (1 - θ)
            end
        end
        if w > 0
            corners[idx] = get(corners, idx, 0.0) + w
        end
    end
    return corners
end

function _grid_point_params(idx::NTuple{N_MODEL_PARAMS,Int})
    return ntuple(i -> PARAM_GRIDS[i][idx[i]], N_MODEL_PARAMS)
end

function _get_or_compute_spectrum(
    energy_sig::UInt64,
    idx::NTuple{N_MODEL_PARAMS,Int},
    energies::AbstractVector{<:Real},
)
    key = (energy_sig, idx)
    grid_params = _grid_point_params(idx)
    cached = lock(CACHE_LOCK) do
        get(SPECTRUM_CACHE, key, nothing)
    end
    if cached !== nothing
        CACHE_HITS[] += 1
        if _verbose_enabled()
            println("GradusXSPEC: using cached spectrum at grid point ($(_format_grid_point(grid_params)))")
        end
        return cached
    end

    CACHE_MISSES[] += 1
    if _verbose_enabled()
        println("GradusXSPEC: evaluating with Gradus at grid point ($(_format_grid_point(grid_params)))")
    end
    spec = _evaluate_model_on_grid(grid_params, energies)
    lock(CACHE_LOCK) do
        return get!(SPECTRUM_CACHE, key, spec)
    end
end

function _interpolated_spectrum(params::NTuple{N_MODEL_PARAMS,Float64}, energies::AbstractVector{<:Real})
    corners = _interpolation_corners(params)
    energy_sig = _energy_signature(energies)
    hits_before = CACHE_HITS[]
    misses_before = CACHE_MISSES[]
    output = zeros(Float64, length(energies) - 1)
    for (idx, weight) in corners
        spec = _get_or_compute_spectrum(energy_sig, idx, energies)
        @inbounds for i in eachindex(output)
            output[i] += weight * spec[i]
        end
    end
    if _verbose_enabled()
        hits = CACHE_HITS[] - hits_before
        misses = CACHE_MISSES[] - misses_before
        n_corners = length(corners)
        cache_size = lock(CACHE_LOCK) do
            length(SPECTRUM_CACHE)
        end
        println(
            "GradusXSPEC: interpolated from $n_corners corner(s) ",
            "($hits cached, $misses evaluated); ",
            "cache holds $cache_size spectrum(s) ",
            "(total hits=$(CACHE_HITS[]), total evaluations=$(CACHE_MISSES[]))",
        )
    end
    return output
end

@ccallable function gradusxspec(energy::Ptr{Cdouble}, Nflux::Cint, parameter::Ptr{Cdouble}, spectrum::Cint, flux::Ptr{Cdouble}, fluxError::Ptr{Cdouble}, init::Ptr{Cchar})::Cint
    try
        _enable_verbose_from_init(init)
        params = unsafe_wrap(Array, parameter, N_MODEL_PARAMS)
        energies = unsafe_wrap(Array, energy, Nflux)
        interp_params = ntuple(i -> Float64(params[i]), N_MODEL_PARAMS)
        if _verbose_enabled()
            println("GradusXSPEC: requested parameters ($(_format_grid_point(interp_params)))")
        end
        flux_array = _interpolated_spectrum(interp_params, energies)

        # copy flux_array back into flux
        for i in 1:length(flux_array)
            unsafe_store!(flux, flux_array[i], i)
        end
        if Nflux > length(flux_array)
            unsafe_store!(flux, 0.0, Nflux)
        end

        return 0
    catch e
        @error "Error in gradusxspec" exception=e
        return 1
    end
end

end
