# External [kerrz](https://git.sr.ht/~fjebaker/kerrz) CLI backend for L(g).
# Pipeline: `emissivity` → FITS → `lineprof` → interpolate onto `g_grid`.

const KERRZ_LOCK = ReentrantLock()
const DEFAULT_KERRZ_PATH = expanduser("~/GitHub/kerrz/zig-out/bin/kerrz")
const DEFAULT_KERRZ_PHOTON_INDEX = 2.0
const DEFAULT_KERRZ_ROUT = 400.0
const DEFAULT_KERRZ_NPHOTONS_LAMP = 3000
const DEFAULT_KERRZ_NPHOTONS_RING = 50_000
const DEFAULT_KERRZ_NG = 1000
# Bump when emissivity/lineprof CLI flags or expected FITS/CSV layout change.
const KERRZ_CACHE_FORMAT = "1"
const KERRZ_VERSION_CACHE = Dict{String, String}()
const KERRZ_VERSION_LOCK = ReentrantLock()

function _env_int(name::AbstractString, default::Int)
    raw = get(ENV, name, "")
    isempty(raw) && return default
    try
        return parse(Int, raw)
    catch
        return default
    end
end

function _env_float_kerrz(name::AbstractString, default::Float64)
    raw = get(ENV, name, "")
    isempty(raw) && return default
    try
        return parse(Float64, raw)
    catch
        return default
    end
end

"""
    kerrz_binary_path() -> String

Resolve the kerrz executable. Override with `GRADUSXSPEC_KERRZ`
(absolute path or a name on `PATH`). Default: `~/GitHub/kerrz/zig-out/bin/kerrz`.
"""
function kerrz_binary_path()
    explicit = get(ENV, "GRADUSXSPEC_KERRZ", "")
    if !isempty(explicit)
        path = expanduser(explicit)
        if isfile(path) && isexecutable(path)
            return path
        end
        found = Sys.which(path)
        found === nothing && throw(ErrorException(
            "GRADUSXSPEC_KERRZ=$(explicit) is not an executable file or PATH entry",
        ))
        return found
    end
    if isfile(DEFAULT_KERRZ_PATH) && isexecutable(DEFAULT_KERRZ_PATH)
        return DEFAULT_KERRZ_PATH
    end
    found = Sys.which("kerrz")
    found === nothing && throw(ErrorException(
        "kerrz binary not found; set GRADUSXSPEC_KERRZ or install at $(DEFAULT_KERRZ_PATH)",
    ))
    return found
end

function _kerrz_nthreads()
    return _env_int("GRADUSXSPEC_KERRZ_NTHREADS", max(1, Threads.nthreads()))
end

function _kerrz_nphotons_lamp()
    return _env_int("GRADUSXSPEC_KERRZ_NPHOTONS_LAMP", DEFAULT_KERRZ_NPHOTONS_LAMP)
end

function _kerrz_nphotons_ring()
    return _env_int("GRADUSXSPEC_KERRZ_NPHOTONS_RING", DEFAULT_KERRZ_NPHOTONS_RING)
end

function _kerrz_rout()
    return _env_float_kerrz("GRADUSXSPEC_KERRZ_ROUT", DEFAULT_KERRZ_ROUT)
end

function _kerrz_em_cache_dir()
    explicit = get(ENV, "GRADUSXSPEC_KERRZ_EM_CACHE_DIR", "")
    isempty(explicit) || return explicit
    return joinpath(first(DEPOT_PATH), "gradusxspec", "kerrz_em")
end

"""
    _kerrz_version_string(bin) -> String

Cached `kerrz --version` output for `bin`. Used in on-disk and L(g) cache keys
so upgrades invalidate stale emissivity FITS / line profiles.
"""
function _kerrz_version_string(bin::AbstractString)
    return lock(KERRZ_VERSION_LOCK) do
        get!(KERRZ_VERSION_CACHE, bin) do
            cmd = Cmd(`$bin --version`; ignorestatus = true)
            buf_out = IOBuffer()
            buf_err = IOBuffer()
            proc = run(pipeline(cmd; stdout = buf_out, stderr = buf_err))
            out = strip(String(take!(buf_out)))
            err = strip(String(take!(buf_err)))
            if success(proc) && !isempty(out)
                return out
            end
            # Still fingerprint something stable so cache keys change if probing fails.
            return "unknown:exit=$(proc.exitcode):out=$(out):err=$(err)"
        end
    end
end

function _kerrz_nphotons_for(variant::Symbol)
    if variant == :kerrz_lamppost
        return _kerrz_nphotons_lamp()
    elseif variant == :kerrz_ring
        return _kerrz_nphotons_ring()
    end
    throw(ArgumentError("unsupported kerrz variant for photon count: $variant"))
end

"""
Fingerprint of kerrz binary + CLI settings that affect L(g) (and emissivity FITS).
"""
function _kerrz_runtime_fingerprint(variant::Symbol; bin::AbstractString = kerrz_binary_path())
    return join(
        [
            abspath(bin),
            _kerrz_version_string(bin),
            KERRZ_CACHE_FORMAT,
            string(variant),
            string(_kerrz_nphotons_for(variant)),
            string(DEFAULT_KERRZ_PHOTON_INDEX),
            string(_kerrz_rout()),
            string(DEFAULT_KERRZ_NG),
            "corotate",
        ],
        "|",
    )
end

function _kerrz_run!(bin::AbstractString, args::Vector{String}; label::AbstractString)
    cmd = Cmd(`$bin $args`; ignorestatus = true)
    buf_out = IOBuffer()
    buf_err = IOBuffer()
    proc = run(pipeline(cmd; stdout = buf_out, stderr = buf_err))
    if !success(proc)
        err = String(take!(buf_err))
        out = String(take!(buf_out))
        throw(ErrorException(
            "kerrz $label failed (exit $(proc.exitcode))\n" *
            "command: kerrz $(join(args, " "))\n" *
            "stderr:\n$(err)\n" *
            "stdout:\n$(out)",
        ))
    end
    return nothing
end

function _emissivity_cache_key(
    bin::AbstractString,
    variant::Symbol,
    spin::Float64,
    height::Float64,
    radius::Union{Float64, Nothing},
    nphotons::Int,
)
    parts = [
        _kerrz_runtime_fingerprint(variant; bin = bin),
        string(spin),
        string(height),
        radius === nothing ? "-" : string(radius),
        string(nphotons),
    ]
    return string(hash(join(parts, "|")); base = 16)
end

function _ensure_emissivity_fits!(
    bin::AbstractString,
    variant::Symbol,
    spin::Float64,
    height::Float64,
    radius::Union{Float64, Nothing},
    dest::AbstractString,
)
    nphotons = _kerrz_nphotons_for(variant)
    cache_dir = _kerrz_em_cache_dir()
    mkpath(cache_dir)
    key = _emissivity_cache_key(bin, variant, spin, height, radius, nphotons)
    cached = joinpath(cache_dir, "em_$(key).fits")
    if isfile(cached) && filesize(cached) > 0
        cp(cached, dest; force = true)
        return dest
    end

    args = String[
        "emissivity",
        "--spin",
        string(spin),
        "--photon-index",
        string(DEFAULT_KERRZ_PHOTON_INDEX),
        "--nphotons",
        string(nphotons),
        "--nthreads",
        string(_kerrz_nthreads()),
        "--output",
        dest,
    ]
    if variant == :kerrz_lamppost
        push!(args, "--lamppost", "h:$(height),vr:0")
    elseif variant == :kerrz_ring
        radius === nothing && throw(ArgumentError("ring corona requires radius"))
        push!(args, "--ring-like", "h:$(height),x:$(radius)")
        push!(args, "--velocity", "corotate")
    else
        throw(ArgumentError("unsupported kerrz emissivity variant: $variant"))
    end

    _kerrz_run!(bin, args; label = "emissivity")
    isfile(dest) || throw(ErrorException("kerrz emissivity did not write $dest"))
    try
        cp(dest, cached; force = true)
    catch
        # Cache write is best-effort.
    end
    return dest
end

function _parse_kerrz_lineprofile(path::AbstractString)
    g = Float64[]
    flux = Float64[]
    open(path, "r") do io
        for (i, line) in enumerate(eachline(io))
            s = strip(line)
            isempty(s) && continue
            if i == 1 && (startswith(lowercase(s), "g") || occursin("flux", lowercase(s)))
                continue
            end
            parts = split(s, ',')
            length(parts) >= 2 || continue
            push!(g, parse(Float64, strip(parts[1])))
            push!(flux, parse(Float64, strip(parts[2])))
        end
    end
    length(g) >= 2 || throw(ErrorException("kerrz lineprof produced too few points in $path"))
    return g, flux
end

function _interp_clamp_nonneg(
    g_src::AbstractVector{<:Real},
    f_src::AbstractVector{<:Real},
    g_dst::AbstractVector{<:Real},
)
    n = length(g_dst)
    out = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        gi = Float64(g_dst[i])
        if gi <= Float64(g_src[1])
            out[i] = max(0.0, Float64(f_src[1]))
        elseif gi >= Float64(g_src[end])
            out[i] = max(0.0, Float64(f_src[end]))
        else
            j = searchsortedlast(g_src, gi)
            j = clamp(j, 1, length(g_src) - 1)
            g0 = Float64(g_src[j])
            g1 = Float64(g_src[j + 1])
            f0 = Float64(f_src[j])
            f1 = Float64(f_src[j + 1])
            t = g1 == g0 ? 0.0 : (gi - g0) / (g1 - g0)
            out[i] = max(0.0, (1 - t) * f0 + t * f1)
        end
    end
    return out
end

function _run_kerrz_lineprof!(
    bin::AbstractString,
    spin::Float64,
    inclination::Float64,
    em_fits::AbstractString,
    line_path::AbstractString,
)
    args = String[
        "lineprof",
        "--spin",
        string(spin),
        "--incl",
        string(inclination),
        "--emissivity-profile",
        em_fits,
        "--rout",
        string(_kerrz_rout()),
        "--ng",
        string(DEFAULT_KERRZ_NG),
        "--nthreads",
        string(_kerrz_nthreads()),
        "--output",
        line_path,
    ]
    _kerrz_run!(bin, args; label = "lineprof")
    isfile(line_path) || throw(ErrorException("kerrz lineprof did not write $line_path"))
    return nothing
end

"""
    kerrz_raw_line_profile(params, g_bins, variant) -> Vector{Float64}

Compute an (approximately) unit-area L(g) via the kerrz CLI and resample onto
`g_bins`. `variant` is `:kerrz_lamppost` (params: spin, inc, h) or
`:kerrz_ring` (params: spin, inc, r, h).
"""
function kerrz_raw_line_profile(
    params::NTuple{N, Float64},
    g_bins::AbstractVector{<:Real},
    variant::Symbol,
) where {N}
    if variant == :kerrz_lamppost
        length(params) == 3 || throw(ArgumentError(
            "kerrz lamppost expects (spin, inc, h); got N=$N",
        ))
        spin, inclination, height = params
        radius = nothing
    elseif variant == :kerrz_ring
        length(params) == 4 || throw(ArgumentError(
            "kerrz ring expects (spin, inc, r, h); got N=$N",
        ))
        spin, inclination, radius, height = params
    else
        throw(ArgumentError("unsupported kerrz variant: $variant"))
    end

    lock(KERRZ_LOCK) do
        bin = kerrz_binary_path()
        mktempdir() do tmp
            em_fits = joinpath(tmp, "em.fits")
            line_path = joinpath(tmp, "line.dat")
            _ensure_emissivity_fits!(bin, variant, spin, height, radius, em_fits)
            _run_kerrz_lineprof!(bin, spin, inclination, em_fits, line_path)
            g_src, f_src = _parse_kerrz_lineprofile(line_path)
            return _interp_clamp_nonneg(g_src, f_src, g_bins)
        end
    end
end
