# Shared approximate memory budget for *bounded* GradusXSPEC caches
# (convolution matrices, line spectra, ring emissivity). Line-profile kernels
# L(g) use a separate unlimited RAM/disk cache and are not counted here.
# Default 16 GiB; override with GRADUSXSPEC_CACHE_LIMIT_GB (set 0 for unlimited).
# The limit is read from the environment at use time so XSPEC sessions pick up
# the value set before `lmod` / process start (not at PackageCompiler build time).

const DEFAULT_CACHE_LIMIT_GIB = 16

function _parse_cache_limit_bytes()
    raw = get(ENV, "GRADUSXSPEC_CACHE_LIMIT_GB", string(DEFAULT_CACHE_LIMIT_GIB))
    gb = try
        parse(Float64, raw)
    catch
        Float64(DEFAULT_CACHE_LIMIT_GIB)
    end
    gb <= 0 && return UInt64(0)
    return UInt64(round(gb * (UInt64(1) << 30)))
end

const CACHE_MEMORY_USED_BYTES = Ref{UInt64}(0)
const CACHE_BUDGET_LOCK = ReentrantLock()

# Evictors registered in preferred order (largest entries first).
# Each callback frees one LRU entry and returns bytes freed (0 if empty).
const CACHE_EVICT_CALLBACKS = Function[]

function register_cache_evictor!(f::Function)
    push!(CACHE_EVICT_CALLBACKS, f)
    return nothing
end

function cache_memory_limit_bytes()
    return _parse_cache_limit_bytes()
end

function cache_memory_used_bytes()
    return CACHE_MEMORY_USED_BYTES[]
end

function cache_memory_limit_gib()
    lim = cache_memory_limit_bytes()
    return lim == 0 ? Inf : lim / (UInt64(1) << 30)
end

function cache_memory_used_gib()
    return CACHE_MEMORY_USED_BYTES[] / (UInt64(1) << 30)
end

function _lru_touch!(order::Vector, key)
    i = findfirst(==(key), order)
    if i !== nothing
        deleteat!(order, i)
    end
    push!(order, key)
    return nothing
end

function _array_nbytes(x::AbstractArray)
    return UInt64(sizeof(x))
end

function _summary_nbytes(x)
    return UInt64(Base.summarysize(x))
end

"""
Free LRU entries across registered caches until `needed` additional bytes fit
under the budget (or caches are empty). Caller must hold `CACHE_BUDGET_LOCK`.
"""
function _evict_until_free!(needed::Integer)
    limit = cache_memory_limit_bytes()
    iszero(limit) && return
    need = UInt64(needed)
    while CACHE_MEMORY_USED_BYTES[] + need > limit
        freed = UInt64(0)
        for f in CACHE_EVICT_CALLBACKS
            freed = f()::UInt64
            freed > 0 && break
        end
        freed == 0 && break
    end
    return nothing
end

function _bounded_cache_lookup!(cache_lock::ReentrantLock, dict::Dict, order::Vector, key)
    return lock(cache_lock) do
        val = get(dict, key, nothing)
        if val !== nothing
            _lru_touch!(order, key)
        end
        return val
    end
end

"""
Insert `value` under `key` if absent, respecting the shared memory budget.
Returns the cached value (existing or newly stored).
`nbytes` is the approximate size of `value`.
"""
function _bounded_cache_put!(
    cache_lock::ReentrantLock,
    dict::Dict,
    order::Vector,
    sizes::Dict,
    key,
    value,
    nbytes::Integer,
)
    nb = UInt64(nbytes)
    return lock(CACHE_BUDGET_LOCK) do
        existing = lock(cache_lock) do
            get(dict, key, nothing)
        end
        existing !== nothing && return existing

        _evict_until_free!(nb)

        return lock(cache_lock) do
            if haskey(dict, key)
                return dict[key]
            end
            dict[key] = value
            sizes[key] = nb
            push!(order, key)
            CACHE_MEMORY_USED_BYTES[] += nb
            return value
        end
    end
end

function _evict_one_from!(
    cache_lock::ReentrantLock,
    dict::Dict,
    order::Vector,
    sizes::Dict,
)
    return lock(cache_lock) do
        isempty(order) && return UInt64(0)
        key = popfirst!(order)
        pop!(dict, key)
        nb = pop!(sizes, key)
        CACHE_MEMORY_USED_BYTES[] -= nb
        return nb::UInt64
    end
end
