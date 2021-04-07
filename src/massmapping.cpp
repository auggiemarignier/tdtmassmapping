#include <iostream>
#include <gsl/gsl_rng.h>

#include "proposals.hpp"

int main()
{
    // setup random number generator
    const gsl_rng_type *T;
    gsl_rng *r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    int total = 10;

    for (int i = 0; i < total; i++)
    {
        double u = gsl_rng_uniform(r);
        std::cout << "i=" << i << "\t u=" << u << "\n";
    }

    GlobalSliceMM global;

    return 0;
}