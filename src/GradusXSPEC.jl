module GradusXSPEC

using Base: @ccallable
using Gradus

struct ComputeResult
    sum::Cdouble
    avg::Cdouble
end

@ccallable function gradusxspec(energy::Ptr{Cdouble}, Nflux::Cint, parameter::Ptr{Cdouble}, spectrum::Cint, flux::Ptr{Cdouble}, fluxError::Ptr{Cdouble}, init::Ptr{Cchar})::Cint
    try
        # do nothing for now!
        println("Called gradusxspec -- trying to do nothing!")
        println("Nflux = ", Nflux)
        # print out entries in the energy array
        energies = unsafe_wrap(Array, energy, Nflux)
        for i in 1:Nflux
            println("Energy[", i, "] = ", energies[i])
        end
        # just as a test put energy^2 in the flux array
        flux_array = unsafe_wrap(Array, flux, Nflux)
        for i in 1:Nflux
            flux_array[i] = energies[i]^2
        end
        return 0
    catch e
        @error "Error in gradusxspec" exception=e
        return 1
    end
end

@ccallable function process_array(ptr::Ptr{Cdouble}, len::Cint, res_ptr::Ptr{ComputeResult})::Cint
    try
        # Wrap the C pointer into a Julia Array (zero-copy)
        # unsafe_wrap(ArrayType, pointer, dimensions)
        data = unsafe_wrap(Array, ptr, len)
        
        # Do "Julia" things (simulated here with basic math)
        total = sum(data)
        average = total / len

	# See if we can call a Gradus function
	m = KerrMetric(M = 1.0, a = 0.998)
	isco = Gradus.isco(m)
	println("ISCO is at ", isco)

        # Store result in the C-provided struct
        # unsafe_store! writes value to the pointer address
        unsafe_store!(res_ptr, ComputeResult(total, average))
        
        return 0 # Return 0 for success
    catch e
        @error "Something went wrong" exception=e
        return 1 # Return 1 for error
    end
end

# Required: Determine real definitions for initialization (PackageCompiler handles most of this)
end
