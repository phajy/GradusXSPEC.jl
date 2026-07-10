module GradusXSPEC

using Base: @ccallable
using Gradus
include("model_definition.jl")
include("line_profile.jl")
include("table_model.jl")
include("convolution.jl")
include("spectrum.jl")

export default_g_grid,
    g_midpoints_from_energy_edges,
    line_profile_kernel,
    line_profile_on_energy_edges,
    E_REST_KEV,
    load_xspec_table,
    get_table,
    interpolate_table_spectrum,
    interpolate_line_kernel,
    build_convolution_matrix,
    convolve_reflection,
    rebin_flux,
    parse_init_string,
    InitConfig,
    evaluate_spectrum,
    evaluate_spectrum_interpolated,
    evaluate_line_spectrum

const LINE_PARAM_GRIDS = ntuple(i -> build_parameter_grid(GRADUS_PARAMETERS[i]), N_GRADUS_PARAMS)
const LINE_CACHE_LOCK = ReentrantLock()
const LINE_SPECTRUM_CACHE = Dict{Tuple{UInt64, NTuple{N_GRADUS_PARAMS, Int}}, Vector{Float64}}()
const VERBOSE = Ref(false)
const LINE_CACHE_HITS = Ref(0)
const LINE_CACHE_MISSES = Ref(0)

function _verbose_enabled()
    return VERBOSE[] || lowercase(get(ENV, "GRADUSXSPEC_VERBOSE", "0")) in ("1", "true", "yes")
end

function _apply_init_config(cfg::InitConfig)
    if cfg.verbose
        VERBOSE[] = true
    end
end

function _apply_init_config(init::Ptr{Cchar})
    _apply_init_config(parse_init_string(init))
end

function _format_line_grid_point(params::NTuple{N_GRADUS_PARAMS, Float64})
    parts = String[]
    for (i, spec) in pairs(GRADUS_PARAMETERS)
        push!(parts, "$(spec.name)=$(params[i])")
    end
    return join(parts, ", ")
end

function _line_grid_bounds_and_weight(value::Float64, grid::Vector{Float64})
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

function _line_interpolation_corners(params::NTuple{N_GRADUS_PARAMS, Float64})
    bounds = ntuple(i -> _line_grid_bounds_and_weight(params[i], LINE_PARAM_GRIDS[i]), N_GRADUS_PARAMS)
    corners = Dict{NTuple{N_GRADUS_PARAMS, Int}, Float64}()
    for mask in 0:(UInt(1) << N_GRADUS_PARAMS) - 1
        idx = ntuple(
            i -> ((mask >> (i - 1)) & UInt(1)) == UInt(1) ? bounds[i][2] : bounds[i][1],
            N_GRADUS_PARAMS,
        )
        weight = 1.0
        for i in 1:N_GRADUS_PARAMS
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

function _line_grid_point_params(idx::NTuple{N_GRADUS_PARAMS, Int})
    return ntuple(i -> LINE_PARAM_GRIDS[i][idx[i]], N_GRADUS_PARAMS)
end

function _get_or_compute_line_spectrum(
    energy_sig::UInt64,
    idx::NTuple{N_GRADUS_PARAMS, Int},
    energies::AbstractVector{<:Real},
)
    key = (energy_sig, idx)
    grid_params = _line_grid_point_params(idx)
    cached = lock(LINE_CACHE_LOCK) do
        get(LINE_SPECTRUM_CACHE, key, nothing)
    end
    if cached !== nothing
        LINE_CACHE_HITS[] += 1
        if _verbose_enabled()
            println("GradusXSPEC: using cached line spectrum at ($(_format_line_grid_point(grid_params)))")
        end
        return cached
    end

    LINE_CACHE_MISSES[] += 1
    if _verbose_enabled()
        println("GradusXSPEC: evaluating line profile at ($(_format_line_grid_point(grid_params)))")
    end
    spec = line_profile_on_energy_edges(energies, grid_params)
    lock(LINE_CACHE_LOCK) do
        return get!(LINE_SPECTRUM_CACHE, key, spec)
    end
end

"""
    evaluate_line_spectrum(energy_edges, gradus_params)

Evaluate the standalone lamppost line profile on `energy_edges` with parameter
grid interpolation and caching.
"""
function evaluate_line_spectrum(
    energy_edges::AbstractVector{<:Real},
    params::NTuple{N_GRADUS_PARAMS, Float64},
)
    corners = _line_interpolation_corners(params)
    energy_sig = _energy_signature(energy_edges)
    output = zeros(Float64, length(energy_edges) - 1)
    for (idx, weight) in corners
        spec = _get_or_compute_line_spectrum(energy_sig, idx, energy_edges)
        @inbounds for i in eachindex(output)
            output[i] += weight * spec[i]
        end
    end
    if _verbose_enabled()
        println(
            "GradusXSPEC: line profile from $(length(corners)) corner(s); ",
            "cache hits=$(LINE_CACHE_HITS[]), misses=$(LINE_CACHE_MISSES[])",
        )
    end
    return output
end

@ccallable function gradusxspec(
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    fluxError::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    try
        cfg = parse_init_string(init)
        _apply_init_config(cfg)
        energies = unsafe_wrap(Array, energy, Nflux)
        params = ntuple(i -> Float64(unsafe_load(parameter, i)), N_XSPEC_FUNC_PARAMS)
        if _verbose_enabled()
            gradus, refl = _split_physics_params(params)
            gradus_text = join(["$(GRADUS_PARAMETERS[i].name)=$(gradus[i])" for i in 1:N_GRADUS_PARAMS], ", ")
            refl_text = join(["$(REFLECTION_PARAMETERS[i].name)=$(refl[i])" for i in 1:N_REFLECTION_PARAMS], ", ")
            println("GradusXSPEC: reflection request ($gradus_text; $refl_text)")
        end
        flux_array = evaluate_spectrum_interpolated(
            energies,
            params;
            table_path = cfg.table_path,
            verbose = _verbose_enabled(),
        )

        for i in 1:length(flux_array)
            unsafe_store!(flux, flux_array[i], i)
        end
        if Nflux > length(flux_array)
            unsafe_store!(flux, 0.0, Nflux)
        end

        return 0
    catch e
        @error "Error in gradusxspec" exception = e
        return 1
    end
end

end
