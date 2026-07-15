const MONITOR_LOCK = ReentrantLock()
# Init-string overrides. `nothing` means "not set via the init string"; the
# getters below then fall back to the environment and finally the defaults.
const MONITOR_PATH = Ref{Union{String, Nothing}}(nothing)
const MONITOR_INTERVAL = Ref{Union{Int, Nothing}}(nothing)
const MONITOR_EVAL_COUNT = Ref(0)
const MONITOR_HIST_BINS = 24
const MONITOR_STATE = Dict{String, Any}()

function _monitor_timestamp()
    return Base.Libc.strftime("%Y-%m-%d %H:%M:%S", round(Int, time()))
end

function _monitor_pad_name(name::AbstractString, width::Int = 12)
    text = string(name)
    return length(text) >= width ? text : text * " " ^ (width - length(text))
end

function _monitor_fmt_num(x::Real)
    return string(x)
end

# GRADUSXSPEC_MONITOR accepts "1"/"true"/"yes" (default path), "0"/"false"/"no"
# (disabled, like verbose parsing), or a custom file path.
function _monitor_env_value()
    env = strip(get(ENV, "GRADUSXSPEC_MONITOR", ""))
    (isempty(env) || lowercase(env) in ("0", "false", "no")) && return nothing
    return env
end

function _monitor_enabled()
    MONITOR_PATH[] !== nothing || _monitor_env_value() !== nothing
end

function _monitor_path()
    if MONITOR_PATH[] !== nothing
        return MONITOR_PATH[]
    end
    env = _monitor_env_value()
    if env === nothing || lowercase(env) in ("1", "true", "yes")
        return DEFAULT_MONITOR_PATH
    end
    return env
end

function _monitor_interval()
    if MONITOR_INTERVAL[] !== nothing
        return MONITOR_INTERVAL[]
    end
    env = strip(get(ENV, "GRADUSXSPEC_MONITOR_INTERVAL", ""))
    if !isempty(env)
        parsed = tryparse(Int, env)
        if parsed !== nothing && parsed >= 1
            return parsed
        end
    end
    return 10
end

# Reset the init-string overrides on every evaluation from the parsed init
# string, mirroring the verbose flag: dropping `monitor` / `monitor_interval`
# from the init string clears the prior setting (the environment variables can
# still enable monitoring via the getters above).
function _apply_monitor_config(cfg::InitConfig)
    MONITOR_PATH[] = cfg.monitor_path
    MONITOR_INTERVAL[] = cfg.monitor_interval
end

function _should_write_monitor_file()
    count = MONITOR_EVAL_COUNT[]
    interval = _monitor_interval()
    return count == 1 || (count % interval == 0)
end

function _monitor_histogram_edges(spec)
    return range(Float64(spec.hard_min), Float64(spec.hard_max), length = MONITOR_HIST_BINS + 1)
end

function _init_monitor_state(rt::ModelRuntime)
    specs = physics_parameters(rt)
    n_params = length(specs)
    edges = ntuple(i -> collect(_monitor_histogram_edges(specs[i])), n_params)
    counts = ntuple(_ -> zeros(Int, MONITOR_HIST_BINS), n_params)
    mins = fill(Inf, n_params)
    maxs = fill(-Inf, n_params)
    return (
        edges = edges,
        counts = counts,
        mins = mins,
        maxs = maxs,
    )
end

function _monitor_state_for(rt::ModelRuntime)
    get!(MONITOR_STATE, rt.definition.name) do
        _init_monitor_state(rt)
    end
end

function _monitor_bin_index(edges::AbstractVector{Float64}, value::Float64)
    idx = searchsortedlast(edges, value)
    return clamp(idx, 1, MONITOR_HIST_BINS)
end

function _record_monitor_params(rt::ModelRuntime, params::Tuple{Vararg{Float64, N}}) where {N}
    state = _monitor_state_for(rt)
    @inbounds for i in 1:N
        v = params[i]
        state.mins[i] = min(state.mins[i], v)
        state.maxs[i] = max(state.maxs[i], v)
        state.counts[i][_monitor_bin_index(state.edges[i], v)] += 1
    end
end

function _ascii_histogram_bar(counts::AbstractVector{Int}; width::Int = 36)
    total = sum(counts)
    isempty(counts) && return "." ^ width
    total == 0 && return "." ^ width
    max_count = maximum(counts)
    chars = Char[]
    for c in counts
        n = max_count == 0 ? 0 : round(Int, width * c / max_count)
        n = c > 0 && n == 0 ? 1 : n
        append!(chars, '#', n)
    end
    bar = String(chars)
    length(bar) >= width && return bar[1:width]
    return bar * "." ^ (width - length(bar))
end

function _format_param_bracket(
    spec,
    grid::Vector{Float64},
    value::Float64,
    width::Int = 32,
)
    lo_idx, hi_idx, θ = _param_grid_bounds_and_weight(value, grid)
    lo = grid[lo_idx]
    hi = grid[hi_idx]
    name = _monitor_pad_name(spec.name)
    if lo_idx == hi_idx
        return "  $(name) $(_monitor_fmt_num(lo))  (on grid)"
    end
    span = hi - lo
    pos = span > 0 ? (value - lo) / span : 0.0
    marker = round(Int, clamp(pos, 0.0, 1.0) * (width - 1)) + 1
    line = fill('-', width)
    line[1] = '['
    line[end] = ']'
    line[marker] = '*'
    θ_text = round(θ, digits = 3)
    return "  $(name) $(_monitor_fmt_num(lo)) .. $(_monitor_fmt_num(hi))  $(String(line))  θ=$(θ_text)"
end

function _format_monitor_histogram(
    spec,
    edges::Vector{Float64},
    counts::Vector{Int},
    hist_lo::Float64,
    hist_hi::Float64,
)
    bar = _ascii_histogram_bar(counts)
    grid_lo = first(edges)
    grid_hi = last(edges)
    total = sum(counts)
    if isfinite(hist_lo)
        history = "history [$(_monitor_fmt_num(hist_lo)) … $(_monitor_fmt_num(hist_hi))]"
    else
        history = "history (none yet)"
    end
    return "  $(spec.name) [$grid_lo … $grid_hi] $history ($total samples)\n  |$bar|"
end

function _write_monitor_file(
    rt::ModelRuntime,
    params::Tuple{Vararg{Float64, N}},
    gradus_corners,
    refl_corners;
    matrix_hits::Int,
    matrix_misses::Int,
    line_hits::Int,
    line_misses::Int,
) where {N}
    lock(MONITOR_LOCK) do
        path = _monitor_path()
        state = _monitor_state_for(rt)
        specs = physics_parameters(rt)
        ng = n_gradus_params(rt)
        gradus_params = ntuple(i -> params[i], ng)
        refl_params = ntuple(i -> params[ng + i], n_reflection_params(rt))

        lines = String[]
        push!(lines, "GradusXSPEC monitor")
        push!(lines, "Model: $(rt.definition.name)")
        push!(lines, "Updated: $(_monitor_timestamp())")
        push!(lines, "Evaluations: $(MONITOR_EVAL_COUNT[]) (file refresh every $(_monitor_interval()) evals)")
        push!(lines, "")

        push!(lines, "Current parameters")
        for (i, spec) in pairs(specs)
            push!(lines, "  $(_monitor_pad_name(spec.name)) $(_monitor_fmt_num(params[i]))")
        end
        push!(lines, "")

        push!(lines, "Gradus interpolation brackets ($(length(gradus_corners)) corner(s))")
        for i in 1:ng
            push!(
                lines,
                _format_param_bracket(
                    rt.definition.gradus_parameters[i],
                    rt.gradus_grids[i],
                    gradus_params[i],
                ),
            )
        end
        push!(lines, "")

        push!(lines, "Reflection interpolation brackets ($(length(refl_corners)) corner(s))")
        for i in 1:n_reflection_params(rt)
            push!(
                lines,
                _format_param_bracket(
                    rt.definition.reflection_parameters[i],
                    rt.reflection_grids[i],
                    refl_params[i],
                ),
            )
        end
        push!(lines, "")

        push!(lines, "Parameter histograms (grid range, fit history, ASCII bar per bin)")
        for i in 1:N
            push!(
                lines,
                _format_monitor_histogram(
                    specs[i],
                    state.edges[i],
                    state.counts[i],
                    state.mins[i],
                    state.maxs[i],
                ),
            )
        end
        push!(lines, "")

        push!(lines, "Cache statistics (cumulative)")
        matrix_total = matrix_hits + matrix_misses
        matrix_rate = matrix_total > 0 ? round(100 * matrix_hits / matrix_total; digits = 1) : 0.0
        push!(
            lines,
            "  Matrix (Gradus blur): hits=$(matrix_hits)  misses=$(matrix_misses)  hit rate=$(matrix_rate)%  cached=$(convolution_matrix_cache_size())",
        )
        line_total = line_hits + line_misses
        line_rate = line_total > 0 ? round(100 * line_hits / line_total; digits = 1) : 0.0
        push!(
            lines,
            "  Line spectrum:      hits=$(line_hits)  misses=$(line_misses)  hit rate=$(line_rate)%  cached=$(line_spectrum_cache_size())",
        )

        open(path, "w") do io
            for line in lines
                println(io, line)
            end
        end
    end
end

function _monitor_after_evaluation(
    rt::ModelRuntime,
    params::Tuple{Vararg{Float64, N}},
    gradus_corners,
    refl_corners,
) where {N}
    _monitor_enabled() || return
    MONITOR_EVAL_COUNT[] += 1
    _record_monitor_params(rt, params)
    _should_write_monitor_file() || return
    matrix_hits, matrix_misses = convolution_matrix_cache_stats()
    line_hits, line_misses = line_spectrum_cache_stats()
    _write_monitor_file(
        rt,
        params,
        gradus_corners,
        refl_corners;
        matrix_hits = matrix_hits,
        matrix_misses = matrix_misses,
        line_hits = line_hits,
        line_misses = line_misses,
    )
end
