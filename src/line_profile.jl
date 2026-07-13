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

function _lamppost_geometry(params::NTuple{4, Float64})
    spin, _, inclination, height = params
    return spin, inclination, height
end

function _lamppost_geometry(params::NTuple{3, Float64})
    spin, inclination, height = params
    return spin, inclination, height
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
    disc_variant::Symbol,
) where {N}
    if disc_variant == :gauss
        length(params) == 1 ||
            throw(ArgumentError("Gaussian kernel expects a single Sigma parameter"))
        return _gaussian_line_profile(params[1], g_bins)
    end

    disc_variant in (:ss, :thin) ||
        throw(ArgumentError("unsupported disc variant: $disc_variant"))
    spin, inclination, height = _lamppost_geometry(params)
    m = KerrMetric(M = 1.0, a = spin)
    d = _build_accretion_disc(m, params, Val(disc_variant))
    x = SVector(0.0, 1000.0, deg2rad(inclination), 0.0)
    corona = LampPostModel(h = height)
    profile = emissivity_profile(m, d, corona; n_samples = 128)
    _, flux = lineprofile(
        m,
        x,
        d,
        profile;
        bins = collect(Float64.(g_bins)),
        method = TransferFunctionMethod(),
        maxrₑ = 400.0,
    )
    return Vector{Float64}(flux)
end

"""
    line_profile_kernel(params, disc_variant; g_grid=default_g_grid()) -> (g, L)

Evaluate a unit-area line-profile kernel `L(g)` for the given lamppost parameters
and accretion-disc variant (`:ss` or `:thin`).
"""
function line_profile_kernel(
    params::NTuple{N, Float64},
    disc_variant::Symbol;
    g_grid::AbstractVector{<:Real} = default_g_grid(),
) where {N}
    g = collect(Float64.(g_grid))
    flux = _raw_line_profile(params, g, disc_variant)
    area = _integrate_g(g, flux)
    area > 0 || error("line profile integral is zero; cannot normalize")
    return g, flux ./ area
end

"""
    line_profile_on_energy_edges(energies, params, disc_variant) -> Vector{Float64}

Evaluate the line profile on the `g` midpoints implied by XSPEC energy edges.
Flux is returned per energy bin and is not forced to unit area.
"""
function line_profile_on_energy_edges(
    energies::AbstractVector{<:Real},
    params::NTuple{N, Float64},
    disc_variant::Symbol,
) where {N}
    g_bins = g_midpoints_from_energy_edges(energies)
    return _raw_line_profile(params, g_bins, disc_variant)
end

# Backwards-compatible overload for the original four-parameter S&S model.
line_profile_kernel(params::NTuple{4, Float64}; kwargs...) =
    line_profile_kernel(params, :ss; kwargs...)
line_profile_on_energy_edges(energies, params::NTuple{4, Float64}) =
    line_profile_on_energy_edges(energies, params, :ss)
