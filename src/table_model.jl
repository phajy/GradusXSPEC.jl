using FITSIO

"""
    XspecTableModel

OGIP/XSPEC additive table model (OGIP 92-009) loaded from a FITS file such as
`xillverD-5.fits`. `SPECTRA` rows follow OGIP ordering (last parameter varies
fastest); `spectra` is reshaped as `(p5, p4, p3, p2, p1, energy)` so parameter
`i` is accessed at index `i` in `view(spectra, i5, i4, i3, i2, i1, :)`.
"""
struct XspecTableModel
    name::String
    param_names::Vector{String}
    param_values::NTuple{5, Vector{Float64}}
    energy_lo::Vector{Float64}
    energy_hi::Vector{Float64}
    spectra::Array{Float64, 6}
end

function _read_parameter_axes(fits::FITS)
    names = String.(strip.(read(fits[2], "NAME")))
    methods = Vector{Int}(read(fits[2], "METHOD"))
    numb = Vector{Int}(read(fits[2], "NUMBVALS"))
    raw = read(fits[2], "VALUE")
    axes = ntuple(5) do i
        Float64.(vec(raw[1:numb[i], i]))
    end
    return names, methods, axes
end

function _resolve_table_path(path::AbstractString)
    if isabspath(path)
        return path
    end
    root = dirname(@__DIR__)
    return joinpath(root, path)
end

"""
    load_xspec_table(path) -> XspecTableModel

Load an XSPEC table model FITS file. Relative paths are resolved from the
GradusXSPEC repository root.
"""
function load_xspec_table(path::AbstractString)
    fits_path = _resolve_table_path(path)
    isfile(fits_path) || error("table model not found: $fits_path")

    f = FITS(fits_path)
    try
        name = String(strip(FITSIO.read_key(f[1], "MODLNAME")[1]))
        param_names, _methods, param_values = _read_parameter_axes(f)
        length(param_names) == 5 || error("expected 5 interpolation parameters, got $(length(param_names))")

        energy_lo = Float64.(read(f[3], "ENERG_LO"))
        energy_hi = Float64.(read(f[3], "ENERG_HI"))
        n_energy = length(energy_lo)

        spec_matrix = transpose(Float64.(read(f[4], "INTPSPEC")))
        n_rows = size(spec_matrix, 1)
        expected_rows = prod(length.(param_values))
        n_rows == expected_rows || error("unexpected SPECTRA row count: $n_rows != $expected_rows")

        n1, n2, n3, n4, n5 = length.(param_values)
        spectra = reshape(spec_matrix, n5, n4, n3, n2, n1, n_energy)

        return XspecTableModel(name, param_names, param_values, energy_lo, energy_hi, spectra)
    finally
        close(f)
    end
end

function energy_midpoints(table::XspecTableModel)
    return (table.energy_lo .+ table.energy_hi) ./ 2
end

function _axis_bounds(values::AbstractVector{<:Real}, x::Float64)
    clamped = clamp(x, first(values), last(values))
    hi = searchsortedfirst(values, clamped)
    if hi <= 1
        return (1, 1, 0.0)
    elseif hi > length(values)
        idx = length(values)
        return (idx, idx, 0.0)
    elseif values[hi] == clamped
        return (hi, hi, 0.0)
    else
        lo = hi - 1
        θ = (clamped - values[lo]) / (values[hi] - values[lo])
        return (lo, hi, θ)
    end
end

function _axis_corners(values::AbstractVector{<:Real}, x::Float64)
    lo, hi, θ = _axis_bounds(values, x)
    if lo == hi
        return Dict(lo => 1.0)
    end
    return Dict(lo => 1 - θ, hi => θ)
end

"""
    interpolate_table_spectrum(table, params) -> Vector{Float64}

Interpolate a rest-frame reflection spectrum at the five table parameter values.
All axes use linear interpolation; `logXi` is already sampled in log-space in the
table.
"""
function interpolate_table_spectrum(table::XspecTableModel, params::NTuple{5, Float64})
    corner_maps = ntuple(i -> _axis_corners(table.param_values[i], params[i]), 5)
    n_energy = size(table.spectra, 6)
    output = zeros(Float64, n_energy)

    for (i1, w1) in corner_maps[1]
        for (i2, w2) in corner_maps[2]
            for (i3, w3) in corner_maps[3]
                for (i4, w4) in corner_maps[4]
                    for (i5, w5) in corner_maps[5]
                        weight = w1 * w2 * w3 * w4 * w5
                        weight == 0 && continue
                        @inbounds output .+= weight .* view(table.spectra, i5, i4, i3, i2, i1, :)
                    end
                end
            end
        end
    end

    return output
end
