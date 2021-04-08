#include <iostream>
#include <gsl/gsl_rng.h>

#include "proposals.hpp"

int main()
{
    int total = 10;

    GlobalSliceMM global("../../../data/Bolshoi_7_clean_256.txt");

    for (int i = 0; i < total; i++)
    {
        double u = global.random.uniform();
        std::cout << "i=" << i << "\t u=" << u << "\n";
    }
    return 0;
}