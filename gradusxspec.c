// XSPEC local-model wrappers for GradusXSPEC.
//
// Build (from repo root, with HEADAS set):
//   ./build-xspec.sh
//
// Load in XSPEC:
//   lmod gradusxspec .
//
// Models defined in model.dat:
//   gradus_lamp_ss   — lamppost + Shakura-Sunyaev thick disc
//   gradus_lamp_thin — lamppost + thin disc
//   test_gauss       — temporary Gaussian blur (narrow ≈ identity)
//
// Note: XSPEC model names may contain underscores, but wrapper function names
// (below) must not.

#include <stdio.h>
#include "build/include/julia_init.h"

extern int graduslampsjxspec(
    const double* energy, int Nflux, const double* parameter, int spectrum,
    double* flux, double* fluxVariance, const char* init);
extern int graduslampthinxspec(
    const double* energy, int Nflux, const double* parameter, int spectrum,
    double* flux, double* fluxVariance, const char* init);
extern int testgaussxspec(
    const double* energy, int Nflux, const double* parameter, int spectrum,
    double* flux, double* fluxVariance, const char* init);

static int julia_initialized = 0;

static void ensure_julia_initialized(void)
{
    if (!julia_initialized) {
        printf("Starting Julia\n");
        init_julia(0, NULL);
        julia_initialized = 1;
    }
}

void graduslampsjulia(
    const double* energy, int Nflux, const double* parameter, int spectrum,
    double* flux, double* fluxVariance, const char* init)
{
    ensure_julia_initialized();
    graduslampsjxspec(energy, Nflux, parameter, spectrum, flux, fluxVariance, init);
}

void graduslampthinjulia(
    const double* energy, int Nflux, const double* parameter, int spectrum,
    double* flux, double* fluxVariance, const char* init)
{
    ensure_julia_initialized();
    graduslampthinxspec(energy, Nflux, parameter, spectrum, flux, fluxVariance, init);
}

void testgaussjulia(
    const double* energy, int Nflux, const double* parameter, int spectrum,
    double* flux, double* fluxVariance, const char* init)
{
    ensure_julia_initialized();
    testgaussxspec(energy, Nflux, parameter, spectrum, flux, fluxVariance, init);
}

// Backwards-compatible wrapper for older model.dat entries.
void gradusjulia(
    const double* energy, int Nflux, const double* parameter, int spectrum,
    double* flux, double* fluxVariance, const char* init)
{
    graduslampsjulia(energy, Nflux, parameter, spectrum, flux, fluxVariance, init);
}
