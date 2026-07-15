"""
    interpolate_line_kernel(g_grid, L, g) -> Float64

Linearly interpolate a line-profile kernel `L(g)` tabulated on `g_grid`.
Values outside the grid are treated as zero (finite kernel support).
"""
function interpolate_line_kernel(
    g_grid::AbstractVector{<:Real},
    L::AbstractVector{<:Real},
    g::Real,
)
    length(g_grid) == length(L) || throw(ArgumentError("g_grid and L must have the same length"))
    isempty(g_grid) && return 0.0

    x = Float64(g)
    x_min = Float64(first(g_grid))
    x_max = Float64(last(g_grid))
    x < x_min && return 0.0
    x > x_max && return 0.0

    hi = searchsortedfirst(g_grid, x)
    if hi <= 1
        return Float64(L[1])
    elseif hi > length(g_grid)
        return Float64(L[end])
    elseif Float64(g_grid[hi]) == x
        return Float64(L[hi])
    else
        lo = hi - 1
        θ = (x - Float64(g_grid[lo])) / (Float64(g_grid[hi]) - Float64(g_grid[lo]))
        return (1 - θ) * Float64(L[lo]) + θ * Float64(L[hi])
    end
end

function _integrate_kernel_over_bins(
    out_lo::Float64,
    out_hi::Float64,
    em_lo::Float64,
    em_hi::Float64,
    g_grid::AbstractVector{<:Real},
    L::AbstractVector{<:Real};
    n_sub::Int = 4,
)
    out_lo > 0 && out_hi > out_lo && em_lo > 0 && em_hi > em_lo || return 0.0
    em_width = em_hi - em_lo
    log_out = range(log(out_lo), log(out_hi), length = n_sub + 1)
    log_em = range(log(em_lo), log(em_hi), length = n_sub + 1)
    integral = 0.0
    @inbounds for io in 1:n_sub
        o_lo = exp(log_out[io])
        o_hi = exp(log_out[io + 1])
        o_mid = sqrt(o_lo * o_hi)
        dO = o_hi - o_lo
        for je in 1:n_sub
            e_lo = exp(log_em[je])
            e_hi = exp(log_em[je + 1])
            e_mid = sqrt(e_lo * e_hi)
            dE = e_hi - e_lo
            g = o_mid / e_mid
            # Photon-number conservation for g = E_obs/E_em with ∫ L(g) dg = 1:
            #   N_obs(E_o) = ∫ L(E_o/E_e) N_em(E_e) dE_e / E_e
            # The 1/E_e is the Jacobian from dg → dE_e (equivalent to dg/g).
            integral += interpolate_line_kernel(g_grid, L, g) * dO * dE / e_mid
        end
    end
    return integral / em_width
end

"""
    build_convolution_matrix(em_lo, em_hi, out_lo, out_hi, g_grid, L; n_sub=4)

Build the matrix `M` for blurring a rest-frame reflection spectrum:

    F = M * R

where `R[j]` is the integrated photon flux in emission bin `j` and `F[i]` is the
integrated photon flux in observed bin `i`. The kernel `L(g)` must have unit
area in `g = E_obs / E_em`. The integrand includes the Jacobian `1/E_em` so that
`∫ L(g) dg = 1` implies photon-number conservation (`∑ F ≈ ∑ R` aside from
photons shifted outside the energy grid).
"""
function build_convolution_matrix(
    em_lo::AbstractVector{<:Real},
    em_hi::AbstractVector{<:Real},
    out_lo::AbstractVector{<:Real},
    out_hi::AbstractVector{<:Real},
    g_grid::AbstractVector{<:Real},
    L::AbstractVector{<:Real};
    n_sub::Int = 4,
)
    length(em_lo) == length(em_hi) || throw(ArgumentError("em_lo and em_hi must match"))
    length(out_lo) == length(out_hi) || throw(ArgumentError("out_lo and out_hi must match"))

    n_out = length(out_lo)
    n_em = length(em_lo)
    M = Matrix{Float64}(undef, n_out, n_em)

    @inbounds for i in 1:n_out
        out_lo_i = Float64(out_lo[i])
        out_hi_i = Float64(out_hi[i])
        for j in 1:n_em
            M[i, j] = _integrate_kernel_over_bins(
                out_lo_i,
                out_hi_i,
                Float64(em_lo[j]),
                Float64(em_hi[j]),
                g_grid,
                L;
                n_sub = n_sub,
            )
        end
    end

    return M
end

"""
    convolve_reflection(R, em_lo, em_hi, g_grid, L; out_lo=em_lo, out_hi=em_hi)

Convolve rest-frame reflection bin fluxes `R` with a unit-area line-profile
kernel `L(g)`.
"""
function convolve_reflection(
    R::AbstractVector{<:Real},
    em_lo::AbstractVector{<:Real},
    em_hi::AbstractVector{<:Real},
    g_grid::AbstractVector{<:Real},
    L::AbstractVector{<:Real};
    out_lo::AbstractVector{<:Real} = em_lo,
    out_hi::AbstractVector{<:Real} = em_hi,
    n_sub::Int = 4,
)
    length(R) == length(em_lo) || throw(ArgumentError("R must match the emission grid"))
    M = build_convolution_matrix(em_lo, em_hi, out_lo, out_hi, g_grid, L; n_sub = n_sub)
    return M * Vector{Float64}(R)
end

"""
    rebin_flux(flux, src_lo, src_hi, dst_lo, dst_hi)

Flux-conserving rebin from source energy bins to destination bins.
"""
function rebin_flux(
    flux::AbstractVector{<:Real},
    src_lo::AbstractVector{<:Real},
    src_hi::AbstractVector{<:Real},
    dst_lo::AbstractVector{<:Real},
    dst_hi::AbstractVector{<:Real},
)
    length(flux) == length(src_lo) || throw(ArgumentError("flux must match source bins"))
    length(dst_lo) == length(dst_hi) || throw(ArgumentError("dst_lo and dst_hi must match"))

    out = zeros(Float64, length(dst_lo))
    @inbounds for k in eachindex(dst_lo)
        dst_lo_k = Float64(dst_lo[k])
        dst_hi_k = Float64(dst_hi[k])
        for j in eachindex(flux)
            src_lo_j = Float64(src_lo[j])
            src_hi_j = Float64(src_hi[j])
            width_j = src_hi_j - src_lo_j
            width_j <= 0 && continue
            overlap = min(src_hi_j, dst_hi_k) - max(src_lo_j, dst_lo_k)
            overlap <= 0 && continue
            out[k] += Float64(flux[j]) * overlap / width_j
        end
    end
    return out
end
