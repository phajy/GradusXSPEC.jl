module GradusXSPEC

using Base: @ccallable
using Gradus
include("model_definition.jl")
include("cache_memory.jl")
include("line_profile.jl")
include("table_model.jl")
include("convolution.jl")
include("blur_grid.jl")
include("spectrum.jl")
include("model_runtime.jl")
include("kernel_cache.jl")
include("monitor.jl")

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
    evaluate_line_spectrum,
    ModelRuntime,
    get_model_runtime,
    build_model_runtimes,
    PACKAGE_NAME,
    ALL_MODELS,
    LAMP_SS_MODEL,
    LAMP_THIN_MODEL,
    RING_THIN_MODEL,
    DISC_THIN_MODEL,
    TEST_GAUSS_MODEL

const LINE_CACHE_LOCK = ReentrantLock()
const LINE_SPECTRUM_CACHE =
    Dict{Tuple{String, UInt64, Tuple{Vararg{Int}}}, Vector{Float64}}()
const LINE_CACHE_ORDER = Any[]
const LINE_CACHE_SIZES = Dict{Any, UInt64}()
const VERBOSE = Ref(false)
const LINE_CACHE_HITS = Ref(0)
const LINE_CACHE_MISSES = Ref(0)

function _evict_one_line_spectrum!()
    return _evict_one_from!(
        LINE_CACHE_LOCK,
        LINE_SPECTRUM_CACHE,
        LINE_CACHE_ORDER,
        LINE_CACHE_SIZES,
    )
end

function _verbose_enabled()
    return VERBOSE[] || lowercase(get(ENV, "GRADUSXSPEC_VERBOSE", "0")) in ("1", "true", "yes")
end

function _apply_init_config(cfg::InitConfig)
    # Track the init string on every call so removing "verbose" from the model
    # init string turns verbose output back off (the environment variable can
    # still force it on via _verbose_enabled).
    VERBOSE[] = cfg.verbose
    _apply_monitor_config(cfg)
end

function _apply_init_config(init::Ptr{Cchar})
    _apply_init_config(parse_init_string(init))
end

function _get_or_compute_line_spectrum(
    rt::ModelRuntime,
    energy_sig::UInt64,
    idx::Tuple{Vararg{Int, N}},
    energies::AbstractVector{<:Real},
) where {N}
    key = (rt.definition.name, energy_sig, idx)
    grid_params = _gradus_grid_point_params(rt, idx)
    cached = _bounded_cache_lookup!(
        LINE_CACHE_LOCK,
        LINE_SPECTRUM_CACHE,
        LINE_CACHE_ORDER,
        key,
    )
    if cached !== nothing
        LINE_CACHE_HITS[] += 1
        if _verbose_enabled()
            println(
                "GradusXSPEC: using cached line spectrum for $(rt.definition.name) at ($(_format_gradus_grid_point(rt, grid_params)))",
            )
        end
        return cached
    end

    LINE_CACHE_MISSES[] += 1
    if _verbose_enabled()
        println(
            "GradusXSPEC: evaluating line profile for $(rt.definition.name) at ($(_format_gradus_grid_point(rt, grid_params)))",
        )
    end
    spec = line_profile_on_energy_edges(
        energies,
        grid_params,
        rt.definition.corona_variant,
        rt.definition.disc_variant,
    )
    return _bounded_cache_put!(
        LINE_CACHE_LOCK,
        LINE_SPECTRUM_CACHE,
        LINE_CACHE_ORDER,
        LINE_CACHE_SIZES,
        key,
        spec,
        _array_nbytes(spec),
    )
end

"""
    evaluate_line_spectrum(energy_edges, gradus_params; model=LAMP_SS_MODEL.name)

Evaluate the standalone lamppost line profile on `energy_edges` with parameter
grid interpolation and caching.
"""
function evaluate_line_spectrum(
    energy_edges::AbstractVector{<:Real},
    params::Tuple{Vararg{Float64, N}};
    model::AbstractString = LAMP_SS_MODEL.name,
) where {N}
    rt = get_model_runtime(model)
    corners = _multilinear_corners(Tuple(rt.gradus_grids), params)
    energy_sig = _energy_signature(energy_edges)
    output = zeros(Float64, length(energy_edges) - 1)
    for (idx, weight) in corners
        spec = _get_or_compute_line_spectrum(rt, energy_sig, idx, energy_edges)
        @inbounds for i in eachindex(output)
            output[i] += weight * spec[i]
        end
    end
    if _verbose_enabled()
        println(
            "GradusXSPEC: $(rt.definition.name) line profile from $(length(corners)) corner(s); ",
            "cache hits=$(LINE_CACHE_HITS[]), misses=$(LINE_CACHE_MISSES[])",
        )
    end
    return output
end

function evaluate_line_spectrum(
    energy_edges::AbstractVector{<:Real},
    params::NTuple{N_GRADUS_PARAMS, Float64},
)
    return evaluate_line_spectrum(energy_edges, params; model = LAMP_SS_MODEL.name)
end

function _xspec_model_entry(
    rt::ModelRuntime,
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    init::Ptr{Cchar},
    flux::Ptr{Cdouble},
)::Cint
    cfg = parse_init_string(init)
    _apply_init_config(cfg)
    # XSPEC passes Nflux energy bins and Nflux+1 bin edges.
    energies = unsafe_wrap(Array, energy, Int(Nflux) + 1)
    n_params = n_physics_params(rt)
    params = ntuple(i -> Float64(unsafe_load(parameter, i)), n_params)
    if _verbose_enabled()
        gradus, refl = _split_physics_params(rt, params)
        gradus_text = join(
            ["$(rt.definition.gradus_parameters[i].name)=$(gradus[i])" for i in 1:n_gradus_params(rt)],
            ", ",
        )
        refl_text = join(
            [
                "$(rt.definition.reflection_parameters[i].name)=$(refl[i])"
                for i in 1:n_reflection_params(rt)
            ],
            ", ",
        )
        println("GradusXSPEC: $(rt.definition.name) reflection request ($gradus_text; $refl_text)")
    end
    flux_array = evaluate_spectrum_interpolated(
        rt,
        energies,
        params;
        table_path = cfg.table_path,
        verbose = _verbose_enabled(),
    )

    n_flux = Int(Nflux)
    length(flux_array) == n_flux || error(
        "GradusXSPEC: expected $n_flux flux bins, got $(length(flux_array))",
    )
    for i in 1:n_flux
        unsafe_store!(flux, flux_array[i], i)
    end

    return 0
end

function _xspec_model_entry_catch(
    rt::ModelRuntime,
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    fluxError::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    try
        return _xspec_model_entry(rt, energy, Nflux, parameter, init, flux)
    catch e
        @error "Error in $(rt.definition.cc_name)" exception = e
        return 1
    end
end

@ccallable function graduslampsjxspec(
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    fluxError::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    return _xspec_model_entry_catch(
        MODEL_RUNTIMES[LAMP_SS_MODEL.name],
        energy,
        Nflux,
        parameter,
        spectrum,
        flux,
        fluxError,
        init,
    )
end

@ccallable function graduslampthinxspec(
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    fluxError::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    return _xspec_model_entry_catch(
        MODEL_RUNTIMES[LAMP_THIN_MODEL.name],
        energy,
        Nflux,
        parameter,
        spectrum,
        flux,
        fluxError,
        init,
    )
end

@ccallable function gradusringthinxspec(
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    fluxError::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    return _xspec_model_entry_catch(
        MODEL_RUNTIMES[RING_THIN_MODEL.name],
        energy,
        Nflux,
        parameter,
        spectrum,
        flux,
        fluxError,
        init,
    )
end

@ccallable function gradusdiscthinxspec(
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    fluxError::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    return _xspec_model_entry_catch(
        MODEL_RUNTIMES[DISC_THIN_MODEL.name],
        energy,
        Nflux,
        parameter,
        spectrum,
        flux,
        fluxError,
        init,
    )
end

@ccallable function testgaussxspec(
    energy::Ptr{Cdouble},
    Nflux::Cint,
    parameter::Ptr{Cdouble},
    spectrum::Cint,
    flux::Ptr{Cdouble},
    fluxError::Ptr{Cdouble},
    init::Ptr{Cchar},
)::Cint
    return _xspec_model_entry_catch(
        MODEL_RUNTIMES[TEST_GAUSS_MODEL.name],
        energy,
        Nflux,
        parameter,
        spectrum,
        flux,
        fluxError,
        init,
    )
end

function line_spectrum_cache_size()
    return length(LINE_SPECTRUM_CACHE)
end

function line_spectrum_cache_stats()
    return LINE_CACHE_HITS[], LINE_CACHE_MISSES[]
end

# Prefer evicting large convolution matrices first, then ring emissivity, then
# line spectra. Shared budget defaults to 16 GiB (GRADUSXSPEC_CACHE_LIMIT_GB).
empty!(CACHE_EVICT_CALLBACKS)
register_cache_evictor!(_evict_one_matrix!)
register_cache_evictor!(_evict_one_ring_emissivity!)
register_cache_evictor!(_evict_one_line_spectrum!)

end
