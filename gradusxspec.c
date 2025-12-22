// XSPEC will use this to build the gradusxspec model
// initpackage gradus model.dat .
// lmod gradus .

#include <stdio.h>
#include "build/include/julia_init.h" 
// Function exported from libGradusXSPEC.dylib (built by PackageCompiler)
extern int gradusxspec(const double* energy, int Nflux, const double* parameter,
    int spectrum, double* flux, double* fluxVariance, const char* init);

static int julia_initialized = 0;

void gradusjulia(const double* energy, int Nflux, const double* parameter,
    int spectrum, double* flux, double* fluxVariance,
    const char* init)
{
    // Initialize Julia only on first call
    if (!julia_initialized) {
        printf("Starting Julia\n");
        init_julia(0, NULL);
        julia_initialized = 1;
    }

    // Initially do absolutely nothing to see if we can get this basic function working!
    printf("Calling Gradus\n");
    gradusxspec(energy, Nflux, parameter, spectrum, flux, fluxVariance, init);

    // Don't shutdown Julia here; keep it running for subsequent calls
    // shutdown_julia(0);
}

// Possble cleanup function to be called at program exit (don't know how to get XSPEC to call this)
// void gradusjulia_cleanup(void)
// {
//     if (julia_initialized) {
//         shutdown_julia(0);
//         julia_initialized = 0;
//     }
// }
