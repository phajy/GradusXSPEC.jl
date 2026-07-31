# Log-energy FFT convolution for multiplicative redshift kernels L(g).
#
# Continuous form (photon density φ):
#   φ_obs(E_o) = ∫ L(E_o/E_e) φ_em(E_e) (dE_e / E_e)
# With τ = ln(E_o/E_e) this is a translation convolution in ln E:
#   φ_obs(e^{u}) = ∫ L(e^τ) φ_em(e^{u-τ}) dτ
# so the FFT convolution theorem applies on a uniform log-E grid.

using FFTW

"""
    convolution_method() -> Symbol

Return `:fft` (default) or `:matrix` from `GRADUSXSPEC_CONVOLVE`.
Accepted aliases: `fourier` → `:fft`; `direct` / `legacy` → `:matrix`.
"""
function convolution_method()
    raw = lowercase(strip(get(ENV, "GRADUSXSPEC_CONVOLVE", "fft")))
    if raw in ("fft", "fourier", "")
        return :fft
    elseif raw in ("matrix", "direct", "legacy")
        return :matrix
    else
        throw(ArgumentError(
            "GRADUSXSPEC_CONVOLVE must be 'fft' or 'matrix' (got $(repr(raw)))",
        ))
    end
end

function _fft_n_log_bins(E_min::Float64, E_max::Float64)
    env = strip(get(ENV, "GRADUSXSPEC_FFT_NBINS", ""))
    if !isempty(env)
        n = tryparse(Int, env)
        n !== nothing && n >= 16 && return n
    end
    # ~0.5% relative spacing; floor at 1024 so narrow kernels match the matrix path.
    n = ceil(Int, log(E_max / E_min) / log(1.005))
    return max(1024, n)
end

"""
Uniform log-energy bin edges spanning `[E_min, E_max]` with `n_bins` bins.
"""
function log_energy_bin_edges(E_min::Float64, E_max::Float64, n_bins::Int)
    E_min > 0 && E_max > E_min || throw(ArgumentError("need 0 < E_min < E_max"))
    n_bins >= 2 || throw(ArgumentError("n_bins must be ≥ 2"))
    edges = exp.(range(log(E_min), log(E_max), length = n_bins + 1))
    return edges[1:end-1], edges[2:end]
end

"""
Energy envelope for FFT blur: span the emission and output grids (the coarse
blur grid already includes g-support pads). Do not extend further with `/g`,
or flux shifted into that extension is dropped when rebinned back to `out`.
"""
function _fft_energy_envelope(
    em_lo::AbstractVector{<:Real},
    em_hi::AbstractVector{<:Real},
    out_lo::AbstractVector{<:Real},
    out_hi::AbstractVector{<:Real},
    g_grid::AbstractVector{<:Real},
)
    E_min = min(Float64(first(em_lo)), Float64(first(out_lo)))
    E_max = max(Float64(last(em_hi)), Float64(last(out_hi)))
    Float64(first(g_grid)) > 0 || throw(ArgumentError("g_grid must be positive"))
    E_min > 0 || throw(ArgumentError("energy grids must be positive"))
    E_max > E_min || throw(ArgumentError("energy envelope is empty"))
    return E_min, E_max
end

"""
Build a length-`n_fft` circular-convolution kernel for spacing `dτ`.

`K[1]` is τ = 0; positive τ occupies the low indices; negative τ wraps to the end.
Values are `L(e^τ) * dτ`.
"""
function _fft_kernel_array(
    g_grid::AbstractVector{<:Real},
    L::AbstractVector{<:Real},
    dτ::Float64,
    n_fft::Int,
)
    K = zeros(Float64, n_fft)
    K[1] = interpolate_line_kernel(g_grid, L, 1.0) * dτ
    half = n_fft ÷ 2
    @inbounds for k in 1:half
        τ = k * dτ
        K[k + 1] = interpolate_line_kernel(g_grid, L, exp(τ)) * dτ
        # Index for -τ: n_fft - k + 1, except k == half on even n when τ = n_fft/2 * dτ
        # is the Nyquist bin (shared); write negative only when distinct.
        if k < n_fft - k
            K[n_fft - k + 1] = interpolate_line_kernel(g_grid, L, exp(-τ)) * dτ
        end
    end
    return K
end

"""
    convolve_reflection_fft(R, em_lo, em_hi, g_grid, L; out_lo, out_hi, n_log)

Convolve rest-frame bin fluxes `R` with `L(g)` via FFT on a uniform log-E grid.
"""
function convolve_reflection_fft(
    R::AbstractVector{<:Real},
    em_lo::AbstractVector{<:Real},
    em_hi::AbstractVector{<:Real},
    g_grid::AbstractVector{<:Real},
    L::AbstractVector{<:Real};
    out_lo::AbstractVector{<:Real} = em_lo,
    out_hi::AbstractVector{<:Real} = em_hi,
    n_log::Union{Nothing, Int} = nothing,
)
    length(R) == length(em_lo) || throw(ArgumentError("R must match the emission grid"))
    length(em_lo) == length(em_hi) || throw(ArgumentError("em_lo and em_hi must match"))
    length(out_lo) == length(out_hi) || throw(ArgumentError("out_lo and out_hi must match"))
    length(g_grid) == length(L) || throw(ArgumentError("g_grid and L must match"))

    E_min, E_max = _fft_energy_envelope(em_lo, em_hi, out_lo, out_hi, g_grid)
    n_bins = n_log === nothing ? _fft_n_log_bins(E_min, E_max) : Int(n_log)
    n_bins = max(n_bins, 16)

    log_lo, log_hi = log_energy_bin_edges(E_min, E_max, n_bins)
    dτ = log(Float64(log_hi[1]) / Float64(log_lo[1]))  # uniform in log E

    R_log = rebin_flux(R, em_lo, em_hi, log_lo, log_hi)
    φ = Vector{Float64}(undef, n_bins)
    @inbounds for i in 1:n_bins
        dE = Float64(log_hi[i]) - Float64(log_lo[i])
        φ[i] = dE > 0 ? Float64(R_log[i]) / dE : 0.0
    end

    # Zero-pad enough for kernel support in τ = ln g.
    g_min = Float64(first(g_grid))
    g_max = Float64(last(g_grid))
    n_kernel = ceil(Int, (log(g_max) - log(g_min)) / dτ) + 3
    n_fft = nextprod([2, 3, 5], n_bins + n_kernel)

    φ_pad = zeros(Float64, n_fft)
    copyto!(φ_pad, 1, φ, 1, n_bins)
    K = _fft_kernel_array(g_grid, L, dτ, n_fft)

    Φ = rfft(φ_pad)
    Κ = rfft(K)
    @inbounds for i in eachindex(Φ)
        Φ[i] *= Κ[i]
    end
    φ_obs_pad = irfft(Φ, n_fft)

    F_log = Vector{Float64}(undef, n_bins)
    @inbounds for i in 1:n_bins
        dE = Float64(log_hi[i]) - Float64(log_lo[i])
        # Numerical noise can be tiny negative.
        F_log[i] = max(φ_obs_pad[i], 0.0) * dE
    end

    return rebin_flux(F_log, log_lo, log_hi, out_lo, out_hi)
end

"""
    convolve_reflection_matrix(R, em_lo, em_hi, g_grid, L; out_lo, out_hi, n_sub)

Legacy bin-integrated matrix convolution (same as the historical `convolve_reflection`).
"""
function convolve_reflection_matrix(
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
