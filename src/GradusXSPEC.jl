module GradusXSPEC

using Base: @ccallable
using Gradus

@ccallable function gradusxspec(energy::Ptr{Cdouble}, Nflux::Cint, parameter::Ptr{Cdouble}, spectrum::Cint, flux::Ptr{Cdouble}, fluxError::Ptr{Cdouble}, init::Ptr{Cchar})::Cint
    try
        params = unsafe_wrap(Array, parameter, 4)
        energies = unsafe_wrap(Array, energy, Nflux)
        flux_array = unsafe_wrap(Array, flux, Nflux)

        println("spin            = ", params[1])
        println("Eddington ratio = ", params[2])
        println("inclination     = ", params[3])
        println("height          = ", params[4])

        # calculate the line profile
        m = KerrMetric(M = 1.0, a = params[1])
        d = ShakuraSunyaev(m, eddington_ratio = params[2])
        θ = params[3]
        x = SVector(0.0, 1000.0, deg2rad(θ), 0.0)
        model = LampPostModel(h = params[4])
        profile = emissivity_profile(m, d, model; n_samples=128)

        e_bins = [(energies[i] + energies[i+1]) / 2 for i in 1:(length(energies)-1)]
        e_bins = e_bins ./ 6.4
        println("Calculating line profile...")
        e_bins, flux_array = lineprofile(m, x, d, profile; bins = e_bins)

        # do we want bin integrated fluxes?
        println("Maximum value of flux_array: ", maximum(flux_array))

        # copy flux_array back into flux
        for i in 1:length(flux_array)
            unsafe_store!(flux, flux_array[i], i)
        end

        return 0
    catch e
        @error "Error in gradusxspec" exception=e
        return 1
    end
end

end
