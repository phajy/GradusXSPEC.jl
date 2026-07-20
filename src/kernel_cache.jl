# Unlimited in-memory (+ optional on-disk) cache of unit-area line-profile
# kernels L(g). Not counted toward GRADUSXSPEC_CACHE_LIMIT_GB.

const KERNEL_CACHE_FORMAT_VERSION = UInt32(1)
const KERNEL_CACHE_MAGIC = b"GXKL"

const KERNEL_CACHE_LOCK = ReentrantLock()
const LINE_KERNEL_CACHE =
    Dict{Tuple{String, UInt64, Tuple{Vararg{Int}}}, Vector{Float64}}()
const KERNEL_CACHE_HITS = Ref(0)
const KERNEL_CACHE_MISSES = Ref(0)
const KERNEL_DISK_LOADS = Ref(0)
const KERNEL_DISK_SAVES = Ref(0)

function _kernel_g_grid_signature(g_grid::AbstractVector{<:Real})
    return UInt64(hash(g_grid))
end

function _kernel_disk_enabled()
    raw = lowercase(get(ENV, "GRADUSXSPEC_KERNEL_CACHE", "1"))
    return raw in ("1", "true", "yes")
end

function _kernel_cache_dir()
    explicit = get(ENV, "GRADUSXSPEC_KERNEL_CACHE_DIR", "")
    isempty(explicit) || return explicit
    return joinpath(first(DEPOT_PATH), "gradusxspec", "kernels")
end

function _kernel_disk_path(key)
    h = hash(key)
    return joinpath(_kernel_cache_dir(), "L_$(string(h; base = 16)).bin")
end

function _g_grids_match(g_cached::AbstractVector, g_grid::AbstractVector)
    length(g_cached) == length(g_grid) || return false
    @inbounds for i in eachindex(g_grid)
        Float64(g_cached[i]) == Float64(g_grid[i]) || return false
    end
    return true
end

function _try_load_kernel_from_disk(key, g_grid::AbstractVector{<:Real})
    _kernel_disk_enabled() || return nothing
    path = _kernel_disk_path(key)
    isfile(path) || return nothing
    try
        open(path, "r") do io
            magic = read(io, 4)
            magic == KERNEL_CACHE_MAGIC || return nothing
            version = read(io, UInt32)
            version == KERNEL_CACHE_FORMAT_VERSION || return nothing
            n = Int(read(io, UInt64))
            n == length(g_grid) || return nothing
            g_cached = Vector{Float64}(undef, n)
            L = Vector{Float64}(undef, n)
            read!(io, g_cached)
            read!(io, L)
            _g_grids_match(g_cached, g_grid) || return nothing
            KERNEL_DISK_LOADS[] += 1
            return L
        end
    catch
        return nothing
    end
end

function _save_kernel_to_disk!(key, g_grid::AbstractVector{<:Real}, L::Vector{Float64})
    _kernel_disk_enabled() || return
    length(L) == length(g_grid) || return
    dir = _kernel_cache_dir()
    try
        mkpath(dir)
        path = _kernel_disk_path(key)
        tmp = path * ".tmp"
        open(tmp, "w") do io
            write(io, KERNEL_CACHE_MAGIC)
            write(io, KERNEL_CACHE_FORMAT_VERSION)
            write(io, UInt64(length(L)))
            write(io, Float64.(g_grid))
            write(io, L)
        end
        mv(tmp, path; force = true)
        KERNEL_DISK_SAVES[] += 1
    catch
        # Disk cache is best-effort; ignore I/O failures.
    end
    return nothing
end

function _kernel_cache_put_ram!(key, L::Vector{Float64})
    lock(KERNEL_CACHE_LOCK) do
        LINE_KERNEL_CACHE[key] = L
    end
    return L
end

function line_kernel_cache_size()
    return length(LINE_KERNEL_CACHE)
end

function line_kernel_cache_stats()
    return KERNEL_CACHE_HITS[], KERNEL_CACHE_MISSES[]
end

function line_kernel_disk_stats()
    return KERNEL_DISK_LOADS[], KERNEL_DISK_SAVES[]
end

"""
    _get_or_compute_line_kernel(rt, gradus_idx; g_grid) -> Vector{Float64}

Unit-area line-profile kernel `L(g)` for a Gradus grid corner. Unlimited RAM
cache; optional disk persistence via `GRADUSXSPEC_KERNEL_CACHE` /
`GRADUSXSPEC_KERNEL_CACHE_DIR`.
"""
function _get_or_compute_line_kernel(
    rt::ModelRuntime,
    gradus_idx::Tuple{Vararg{Int, N}};
    g_grid::AbstractVector{<:Real} = default_g_grid(),
) where {N}
    g_sig = _kernel_g_grid_signature(g_grid)
    key = (rt.definition.name, g_sig, gradus_idx)

    cached = lock(KERNEL_CACHE_LOCK) do
        get(LINE_KERNEL_CACHE, key, nothing)
    end
    if cached !== nothing
        KERNEL_CACHE_HITS[] += 1
        return cached
    end

    from_disk = _try_load_kernel_from_disk(key, g_grid)
    if from_disk !== nothing
        return _kernel_cache_put_ram!(key, from_disk)
    end

    KERNEL_CACHE_MISSES[] += 1
    gradus_params = _gradus_grid_point_params(rt, gradus_idx)
    _, L = line_profile_kernel(
        gradus_params,
        rt.definition.corona_variant,
        rt.definition.disc_variant;
        g_grid = g_grid,
    )
    L64 = Vector{Float64}(L)
    _kernel_cache_put_ram!(key, L64)
    _save_kernel_to_disk!(key, g_grid, L64)
    return L64
end
