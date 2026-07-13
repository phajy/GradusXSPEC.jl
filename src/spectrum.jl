"""
    InitConfig

Parsed XSPEC model init string. The first token is the reflection table path;
`verbose` may appear as an additional token. `monitor` or `monitor=/path/to/file`
enables fit diagnostics written to a file; `monitor_interval=N` controls how
often the file is refreshed (default 10 evaluations).
"""
struct InitConfig
    table_path::String
    verbose::Bool
    monitor_path::Union{String, Nothing}
    monitor_interval::Union{Int, Nothing}
end

function _parse_monitor_token(token::AbstractString)
    key, _, value = rpartition(token, "=")
    if lowercase(key) == "monitor"
        path = strip(value)
        return isempty(path) ? DEFAULT_MONITOR_PATH : path
    end
    return nothing
end

function _parse_monitor_interval_token(token::AbstractString)
    key, _, value = rpartition(token, "=")
    if lowercase(key) == "monitor_interval"
        parsed = tryparse(Int, strip(value))
        (parsed === nothing || parsed < 1) && return nothing
        return parsed
    end
    return nothing
end

function parse_init_string(init::AbstractString)
    tokens = filter(!isempty, split(strip(init)))
    table_path = DEFAULT_TABLE_PATH
    verbose = false
    monitor_path = nothing
    monitor_interval = nothing
    for token in tokens
        lowered = lowercase(token)
        if lowered == "verbose"
            verbose = true
        elseif startswith(lowered, "monitor_interval")
            monitor_interval = _parse_monitor_interval_token(token)
        elseif startswith(lowered, "monitor")
            monitor_path = _parse_monitor_token(token)
        elseif token != "0"
            table_path = token
        end
    end
    return InitConfig(table_path, verbose, monitor_path, monitor_interval)
end

function parse_init_string(init::Ptr{Cchar})
    init == C_NULL && return InitConfig(DEFAULT_TABLE_PATH, false, nothing, nothing)
    return parse_init_string(strip(unsafe_string(init)))
end

const TABLE_LOCK = ReentrantLock()
const LOADED_TABLE = Ref{Union{XspecTableModel,Nothing}}(nothing)
const LOADED_TABLE_PATH = Ref("")

function get_table(table_path::AbstractString = DEFAULT_TABLE_PATH)
    resolved = _resolve_table_path(table_path)
    lock(TABLE_LOCK) do
        if LOADED_TABLE[] === nothing || LOADED_TABLE_PATH[] != resolved
            LOADED_TABLE[] = load_xspec_table(resolved)
            LOADED_TABLE_PATH[] = resolved
        end
        return LOADED_TABLE[]
    end
end

function _energy_signature(energies::AbstractVector{<:Real})
    h = UInt(0)
    @inbounds for e in energies
        h = hash(Float64(e), h)
    end
    return UInt64(h)
end

function _rebin_to_energy_edges(
    flux::AbstractVector{<:Real},
    src_lo::AbstractVector{<:Real},
    src_hi::AbstractVector{<:Real},
    energy_edges::AbstractVector{<:Real},
)
    length(energy_edges) >= 2 || throw(ArgumentError("energy_edges must contain at least two values"))
    dst_lo = Float64.(energy_edges[1:(end - 1)])
    dst_hi = Float64.(energy_edges[2:end])
    return rebin_flux(flux, src_lo, src_hi, dst_lo, dst_hi)
end

function _param_grid_bounds_and_weight(value::Float64, grid::Vector{Float64})
    clamped = clamp(value, first(grid), last(grid))
    hi = searchsortedfirst(grid, clamped)
    if hi <= 1
        return (1, 1, 0.0)
    elseif hi > length(grid)
        idx = length(grid)
        return (idx, idx, 0.0)
    elseif grid[hi] == clamped
        return (hi, hi, 0.0)
    else
        lo = hi - 1
        θ = (clamped - grid[lo]) / (grid[hi] - grid[lo])
        return (lo, hi, θ)
    end
end

function _multilinear_corners(
    param_grids::NTuple{N, Vector{Float64}},
    params::NTuple{N, Float64},
) where {N}
    bounds = ntuple(i -> _param_grid_bounds_and_weight(params[i], param_grids[i]), N)
    corners = Dict{NTuple{N, Int}, Float64}()
    for mask in 0:(UInt(1) << N) - 1
        idx = ntuple(
            i -> ((mask >> (i - 1)) & UInt(1)) == UInt(1) ? bounds[i][2] : bounds[i][1],
            N,
        )
        weight = 1.0
        for i in 1:N
            θ = bounds[i][3]
            if ((mask >> (i - 1)) & UInt(1)) == UInt(1)
                weight *= θ
            else
                weight *= (1 - θ)
            end
        end
        if weight > 0
            corners[idx] = get(corners, idx, 0.0) + weight
        end
    end
    return corners
end
