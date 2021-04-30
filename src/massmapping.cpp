#include <iostream>
#include <gsl/gsl_rng.h>

#include "birthslice.hpp"
#include "deathslice.hpp"
#include "valueslice.hpp"

#include "proposals.hpp"

int main()
{
    int total = 10;

    GlobalSliceMM global("../../../data/Bolshoi_7_clean_256.txt");
    BirthSlice birth(global);
    DeathSlice death(global);
    ValueSlice value(global);

    global.current_likelihood = global.likelihood(global.current_log_normalization);
    printf("Initial Likelihood: %f\n", global.current_likelihood);

    for (int i = 0; i < total; i++)
    {
        double u = global.random.uniform();
        std::cout << "i=" << i << "\t u=" << u << "\n";
    }
    return 0;
}