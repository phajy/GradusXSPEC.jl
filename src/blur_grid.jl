# Working energy grid for relativistic blur / convolution matrices.
#
# Default: core band [2, 150] keV with ΔE = max(0.1 keV, 1% of E), plus
# extremely coarse pads below/above so photons can redshift into/out of the
# core (approximate flux conservation given L(g) support). Override via ENV.

const DEFAULT_BLUR_EMIN_KEV = 2.0
const DEFAULT_BLUR_EMAX_KEV = 150.0
const DEFAULT_BLUR_DE_ABS_KEV = 0.1
const DEFAULT_BLUR_DE_REL = 0.01
# High-energy pad: very wide bins (50% relative or 20 keV, whichever larger).
const DEFAULT_BLUR_PAD_DE_ABS_KEV = 20.0
const DEFAULT_BLUR_PAD_DE_REL = 0.5

function _env_float(name::AbstractString, default::Float64)
    raw = get(ENV, name, "")
    isempty(raw) && return default
    try
        return parse(Float64, raw)
    catch
        return default
    end
end

function _blur_use_native_table_grid()
    raw = lowercase(get(ENV, "GRADUSXSPEC_BLUR_NATIVE", "0"))
    return raw in ("1", "true", "yes")
end

function blur_emin_kev()
    return _env_float("GRADUSXSPEC_BLUR_EMIN", DEFAULT_BLUR_EMIN_KEV)
end

function blur_emax_kev()
    return _env_float("GRADUSXSPEC_BLUR_EMAX", DEFAULT_BLUR_EMAX_KEV)
end

function blur_de_abs_kev()
    return _env_float("GRADUSXSPEC_BLUR_DE_ABS", DEFAULT_BLUR_DE_ABS_KEV)
end

function blur_de_rel()
    return _env_float("GRADUSXSPEC_BLUR_DE_REL", DEFAULT_BLUR_DE_REL)
end

"""
Append core-band edges from `E` up to `Emax` with
`ΔE = max(de_abs, de_rel * E)`. `edges` must already contain the starting `E`.
"""
function _append_energy_edges!(
    edges::Vector{Float64},
    Emax::Float64;
    de_abs::Float64,
    de_rel::Float64,
)
    E = edges[end]
    while E < Emax - 1e-12
        dE = max(de_abs, de_rel * E)
        Enext = min(E + dE, Emax)
        if Enext <= E
            push!(edges, Emax)
            break
        end
        push!(edges, Enext)
        E = Enext
    end
    return edges
end

"""
    build_blur_energy_edges(table_lo, table_hi; g_grid=default_g_grid()) -> edges

Energy edges for the blur working grid. Core `[GRADUSXSPEC_BLUR_EMIN,
GRADUSXSPEC_BLUR_EMAX]` uses fine coarse-binning; pads extend toward
`Emin/g_max` and `Emax/g_min` (clamped to the table range) with extremely
coarse bins. Set `GRADUSXSPEC_BLUR_NATIVE=1` to use the native table grid.
"""
function build_blur_energy_edges(
    table_lo::AbstractVector{<:Real},
    table_hi::AbstractVector{<:Real};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
)
    table_emin = Float64(first(table_lo))
    table_emax = Float64(last(table_hi))

    if _blur_use_native_table_grid()
        edges = Vector{Float64}(undef, length(table_lo) + 1)
        edges[1:end-1] .= Float64.(table_lo)
        edges[end] = table_emax
        return edges
    end

    core_lo = blur_emin_kev()
    core_hi = blur_emax_kev()
    core_lo < core_hi || throw(ArgumentError(
        "blur Emin ($core_lo) must be < Emax ($core_hi)",
    ))

    g_min = Float64(first(g_grid))
    g_max = Float64(last(g_grid))
    g_min > 0 || throw(ArgumentError("g_grid must be positive"))

    # Emission that can land in the core under g ∈ [g_min, g_max].
    pad_lo = max(table_emin, core_lo / g_max)
    pad_hi = min(table_emax, core_hi / g_min)

    # Clamp core to the table / pad envelope.
    core_lo = clamp(core_lo, pad_lo, pad_hi)
    core_hi = clamp(core_hi, core_lo, pad_hi)

    edges = Float64[pad_lo]

    # Low pad: a single wide bin into the core (extremely coarse).
    if pad_lo < core_lo - 1e-12
        push!(edges, core_lo)
    end

    # Core band.
    _append_energy_edges!(
        edges,
        core_hi;
        de_abs = blur_de_abs_kev(),
        de_rel = blur_de_rel(),
    )

    # High pad: extremely coarse bins up to pad_hi.
    if pad_hi > core_hi + 1e-12
        _append_energy_edges!(
            edges,
            pad_hi;
            de_abs = DEFAULT_BLUR_PAD_DE_ABS_KEV,
            de_rel = DEFAULT_BLUR_PAD_DE_REL,
        )
    end

    length(edges) >= 2 || throw(ArgumentError("blur energy grid is empty"))
    return edges
end

function blur_energy_bin_edges(
    table_lo::AbstractVector{<:Real},
    table_hi::AbstractVector{<:Real};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
)
    edges = build_blur_energy_edges(table_lo, table_hi; g_grid = g_grid)
    lo = edges[1:end-1]
    hi = edges[2:end]
    return lo, hi
end

function blur_grid_signature(
    blur_lo::AbstractVector{<:Real},
    blur_hi::AbstractVector{<:Real},
)
    return UInt64(hash(blur_hi, hash(blur_lo)))
end
