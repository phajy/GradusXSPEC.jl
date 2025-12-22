#include <stdio.h>

// 1. Include the standard initialization header provided by PackageCompiler
// Adjust path to match where it actually is (inside build/include)
#include "build/include/julia_init.h" 

// 2. Manually declare your Julia function here.
// PackageCompiler doesn't know if your function takes ints or doubles,
// so you must write this C signature yourself to match the Julia @ccallable one.
int process_array(double *ptr, int len, void *res_ptr);
int gradusxspec(const double* energy, int Nflux, const double* parameter,
    int spectrum, double* flux, double* fluxVariance,
    const char* init);

// Define your struct
typedef struct {
    double sum;
    double avg;
} ComputeResult;

int main() {
    // Initialize Julia (function comes from julia_init.h)
    init_julia(0, NULL);

    double data[] = {10.0, 20.0, 30.0, 40.0};
    int len = 4;
    ComputeResult res;

    printf("Calling Julia...\n");
    
    // Call your function
    int status = process_array(data, len, &res);

    if (status == 0) {
        printf("Success! Sum: %.2f, Avg: %.2f\n", res.sum, res.avg);
    } else {
        printf("Julia reported an error.\n");
    }

    // Call our XSPEC function to test it
    const double energy[] = {1.0, 2.0, 3.0};
    const double parameter[] = {0.5, 1.5};
    double flux[3];
    double fluxVariance[3];
    gradusxspec(energy, 3, parameter, 0, flux, fluxVariance, "init string");
    // Print flux to verify
    for (int i = 0; i < 3; i++) {
        printf("Flux[%d]: %.2f\n", i, flux[i]);
    }

    shutdown_julia(0);
    return 0;
}
