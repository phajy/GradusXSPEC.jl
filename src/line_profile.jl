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

function _raw_line_profile(params::NTuple{4, Float64}, g_bins::AbstractVector{<:Real})
    spin, eddington, inclination, height = params
    m = KerrMetric(M = 1.0, a = spin)
    d = ShakuraSunyaev(m, eddington_ratio = eddington)
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
    line_profile_kernel(params; g_grid=default_g_grid()) -> (g, L)

Evaluate a unit-area line-profile kernel `L(g)` for the four Gradus parameters
`(spin, Eddington, inc, h)`.
"""
function line_profile_kernel(
    params::NTuple{4, Float64};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
)
    g = collect(Float64.(g_grid))
    flux = _raw_line_profile(params, g)
    area = _integrate_g(g, flux)
    area > 0 || error("line profile integral is zero; cannot normalize")
    return g, flux ./ area
end

"""
    line_profile_on_energy_edges(energies, params) -> Vector{Float64}

Evaluate the line profile on the `g` midpoints implied by XSPEC energy edges.
Flux is returned per energy bin and is not forced to unit area (matching the
current standalone XSPEC model).
"""
function line_profile_on_energy_edges(
    energies::AbstractVector{<:Real},
    params::NTuple{4, Float64},
)
    g_bins = g_midpoints_from_energy_edges(energies)
    return _raw_line_profile(params, g_bins)
end
