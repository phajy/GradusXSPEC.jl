using Gradus

const E_REST_KEV = 6.4

"""
    default_g_grid(; n=256, g_min=0.1, g_max=1.5)

Uniformly spaced redshift grid `g = ν_obs / ν_em` for line-profile kernels.
"""
function default_g_grid(; n::Int = 256, g_min::Real = 0.1, g_max::Real = 1.5)
    return collect(range(Float64(g_min), Float64(g_max), length = n))
end

"""
    g_midpoints_from_energy_edges(energies; E_rest=E_REST_KEV)

Convert XSPEC energy edges (keV) to `g` at bin midpoints, assuming a rest-frame
line energy `E_rest`.
"""
function g_midpoints_from_energy_edges(
    energies::AbstractVector{<:Real};
    E_rest::Real = E_REST_KEV,
)
    n_bins = length(energies) - 1
    g = Vector{Float64}(undef, n_bins)
    inv_E_rest = 1.0 / Float64(E_rest)
    @inbounds for i in 1:n_bins
        g[i] = (Float64(energies[i]) + Float64(energies[i + 1])) * 0.5 * inv_E_rest
    end
    return g
end

function _integrate_g(g::AbstractVector{<:Real}, f::AbstractVector{<:Real})
    length(g) == length(f) || throw(ArgumentError("g and f must have the same length"))
    length(g) >= 2 || return 0.0
    area = 0.0
    @inbounds for i in 1:(length(g) - 1)
        dg = g[i + 1] - g[i]
        area += 0.5 * (f[i] + f[i + 1]) * dg
    end
    return area
end

function _build_accretion_disc(m::KerrMetric, params::NTuple{4, Float64}, ::Val{:ss})
    _, eddington, _, _ = params
    return ShakuraSunyaev(m, eddington_ratio = eddington)
end

function _build_accretion_disc(m::KerrMetric, ::NTuple{3, Float64}, ::Val{:thin})
    return ThinDisc(0.0, Inf)
end

function _build_accretion_disc(m::KerrMetric, ::NTuple{4, Float64}, ::Val{:thin})
    return ThinDisc(0.0, Inf)
end

function _lamppost_geometry(params::NTuple{4, Float64})
    spin, _, inclination, height = params
    return spin, inclination, height
end

function _lamppost_geometry(params::NTuple{3, Float64})
    spin, inclination, height = params
    return spin, inclination, height
end

function _ring_geometry(params::NTuple{4, Float64})
    spin, inclination, radius, height = params
    return spin, inclination, radius, height
end

function _build_corona(params::NTuple{N, Float64}, ::Val{:lamppost}) where {N}
    _, _, height = _lamppost_geometry(params)
    return LampPostModel(h = height)
end

function _build_corona(params::NTuple{4, Float64}, ::Val{:ring})
    _, _, radius, height = _ring_geometry(params)
    return RingCorona(; r = radius, h = height)
end

function _corona_n_samples(::Val{:lamppost})
    return 128
end

function _corona_n_samples(::Val{:ring})
    # RingCorona needs n_samples > 2 * extrema_iter (default 80).
    return 256
end

function _corona_n_samples(::Val{:disc})
    # Each stacked RingCorona needs the same sample count as the ring model.
    return 256
end

# Match Gradus DiscCorona defaults: n concentric rings with r·Δr weighting
# (uniform corona surface brightness). Gradus DiscCorona.emissivity_profile is
# currently unusable here (DiscCoronaProfile expects RingCoronaProfile, but the
# optimized RingCorona path returns RingApproximation).
const DISC_CORONA_N_RINGS = 10
const DISC_CORONA_R_INNER = 1e-2


"""
Clamp non-positive branch emissivities in a Gradus `RingApproximation`.

Gradus `emissivity_at` does `log2.(br.ε)` and DomainErrors on slightly negative
numerical ε (see `gradus_bugs.md` / `reproduce_disc_log2.jl`).
"""
function _clamp_ring_approximation_emissivity!(profile)
    profile isa Gradus.RingApproximation || return profile
    for branch_group in profile.branches
        for br in branch_group
            @inbounds for i in eachindex(br.ε)
                if !(br.ε[i] > 0)
                    br.ε[i] = eps(typeof(br.ε[i]))
                end
            end
        end
    end
    return profile
end

"""
Build a disc-corona emissivity function matching Gradus `DiscCorona`:

    radii = range(r_inner, r_outer, n_rings)
    ε(ρ) = Σᵢ ε_ringᵢ(ρ) · rᵢ · Δr
"""
function _disc_corona_emissivity(
    m::KerrMetric,
    d,
    r_outer::Float64,
    height::Float64;
    n_rings::Int = DISC_CORONA_N_RINGS,
    n_samples::Int = _corona_n_samples(Val(:disc)),
    r_inner::Float64 = DISC_CORONA_R_INNER,
)
    r_outer >= r_inner || throw(ArgumentError(
        "disc corona outer radius ($r_outer) must be >= inner ($r_inner)",
    ))
    radii = collect(range(r_inner, r_outer; length = n_rings))
    δr = n_rings > 1 ? (radii[2] - radii[1]) : r_outer
    profiles = map(radii) do r
        profile = emissivity_profile(
            m,
            d,
            RingCorona(; r = r, h = height);
            n_samples = n_samples,
        )
        _clamp_ring_approximation_emissivity!(profile)
    end
    return ρ -> begin
        total = 0.0
        @inbounds for i in eachindex(radii)
            ε = emissivity_at(profiles[i], ρ)
            ε_s = ε isa Number ? Float64(ε) : Float64(first(ε))
            total += ε_s * radii[i] * δr
        end
        return total
    end
end

"""
Scalar emissivity for Gradus `lineprofile`.

`RingApproximation` (from `RingCorona`) returns a length-1 vector from
`emissivity_at(profile, r)` even for scalar `r`, which breaks transfer-function
integration. Lamppost `RadialDiscProfile` already returns a scalar.
"""
function _scalar_emissivity(profile)
    return r -> begin
        ε = emissivity_at(profile, r)
        return ε isa Number ? Float64(ε) : Float64(first(ε))
    end
end

function _gaussian_line_profile(σ::Float64, g_bins::AbstractVector{<:Real})
    g = collect(Float64.(g_bins))
    n = length(g)
    L = zeros(Float64, n)
    n < 1 && return L

    # Narrow Gaussians are handled as a true identity (no convolution) upstream.
    # Sampling a discrete delta here and then integrating it continuously does not
    # preserve flux, so refuse to build a broken kernel for σ ≪ Δg.
    if n >= 2
        dg = (g[end] - g[1]) / (n - 1)
        if σ < 0.5 * dg
            throw(ArgumentError(
                "Gaussian Sigma=$(σ) is narrower than the g-grid spacing ($(dg)); " *
                "use the identity (no-blur) path instead of sampling a discrete delta",
            ))
        end
    end

    σ_eff = max(σ, 1e-12)
    inv_norm = 1.0 / (σ_eff * sqrt(2π))
    @inbounds for i in 1:n
        x = (g[i] - 1.0) / σ_eff
        L[i] = inv_norm * exp(-0.5 * x * x)
    end
    return L
end

"""
    gaussian_is_identity(σ; g_grid=default_g_grid()) -> Bool

True when `σ` is too narrow to resolve on `g_grid`. In that regime the correct
convolution is the identity: return the rebinned table spectrum with no blur.
"""
function gaussian_is_identity(
    σ::Float64;
    g_grid::AbstractVector{<:Real} = default_g_grid(),
)
    length(g_grid) < 2 && return true
    dg = (Float64(g_grid[end]) - Float64(g_grid[1])) / (length(g_grid) - 1)
    return σ < 0.5 * dg
end

function _raw_line_profile(
    params::NTuple{N, Float64},
    g_bins::AbstractVector{<:Real},
    corona_variant::Symbol,
    disc_variant::Symbol,
) where {N}
    if corona_variant == :gauss || disc_variant == :gauss
        length(params) == 1 ||
            throw(ArgumentError("Gaussian kernel expects a single Sigma parameter"))
        return _gaussian_line_profile(params[1], g_bins)
    end

    corona_variant in (:lamppost, :ring, :disc) ||
        throw(ArgumentError("unsupported corona variant: $corona_variant"))
    disc_variant in (:ss, :thin) ||
        throw(ArgumentError("unsupported disc variant: $disc_variant"))

    if corona_variant == :lamppost
        spin, inclination, _ = _lamppost_geometry(params)
    else
        spin, inclination, _, _ = _ring_geometry(params)
    end

    m = KerrMetric(M = 1.0, a = spin)
    d = _build_accretion_disc(m, params, Val(disc_variant))
    x = SVector(0.0, 1000.0, deg2rad(inclination), 0.0)
    bins = collect(Float64.(g_bins))

    ε = if corona_variant == :disc
        _, _, r_outer, height = _ring_geometry(params)
        _disc_corona_emissivity(m, d, r_outer, height)
    else
        corona = _build_corona(params, Val(corona_variant))
        profile = emissivity_profile(
            m,
            d,
            corona;
            n_samples = _corona_n_samples(Val(corona_variant)),
        )
        _clamp_ring_approximation_emissivity!(profile)
        _scalar_emissivity(profile)
    end

    _, flux = lineprofile(
        bins,
        ε,
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        maxrₑ = 400.0,
    )
    return Vector{Float64}(flux)
end

"""
    line_profile_kernel(params, corona_variant, disc_variant; g_grid=default_g_grid()) -> (g, L)

Evaluate a unit-area line-profile kernel `L(g)` for the given corona/disc geometry.
`corona_variant` is `:lamppost`, `:ring`, `:disc`, or `:gauss`; `disc_variant`
is `:ss`, `:thin`, or `:gauss`.
"""
function line_profile_kernel(
    params::NTuple{N, Float64},
    corona_variant::Symbol,
    disc_variant::Symbol;
    g_grid::AbstractVector{<:Real} = default_g_grid(),
) where {N}
    g = collect(Float64.(g_grid))
    flux = _raw_line_profile(params, g, corona_variant, disc_variant)
    area = _integrate_g(g, flux)
    area > 0 || error("line profile integral is zero; cannot normalize")
    return g, flux ./ area
end

# Backwards-compatible: disc_variant alone implies lamppost (or gauss for :gauss).
function line_profile_kernel(
    params::NTuple{N, Float64},
    disc_variant::Symbol;
    kwargs...,
) where {N}
    corona = disc_variant == :gauss ? :gauss : :lamppost
    return line_profile_kernel(params, corona, disc_variant; kwargs...)
end

"""
    line_profile_on_energy_edges(energies, params, corona_variant, disc_variant) -> Vector{Float64}

Evaluate the line profile on the `g` midpoints implied by XSPEC energy edges.
Flux is returned per energy bin and is not forced to unit area.
"""
function line_profile_on_energy_edges(
    energies::AbstractVector{<:Real},
    params::NTuple{N, Float64},
    corona_variant::Symbol,
    disc_variant::Symbol,
) where {N}
    g_bins = g_midpoints_from_energy_edges(energies)
    return _raw_line_profile(params, g_bins, corona_variant, disc_variant)
end

function line_profile_on_energy_edges(
    energies::AbstractVector{<:Real},
    params::NTuple{N, Float64},
    disc_variant::Symbol,
) where {N}
    corona = disc_variant == :gauss ? :gauss : :lamppost
    return line_profile_on_energy_edges(energies, params, corona, disc_variant)
end

# Backwards-compatible overload for the original four-parameter S&S model.
line_profile_kernel(params::NTuple{4, Float64}; kwargs...) =
    line_profile_kernel(params, :lamppost, :ss; kwargs...)
line_profile_on_energy_edges(energies, params::NTuple{4, Float64}) =
    line_profile_on_energy_edges(energies, params, :lamppost, :ss)
