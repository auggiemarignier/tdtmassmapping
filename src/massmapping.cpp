#include <iostream>
#include <gsl/gsl_rng.h>

#include "proposals.hpp"

int main()
{
    int total = 10;
    int seed = 1;

    int kmax = 100;

    double Pb = 0.05;

    int wavelet_xy = 4;
    int degreex = 8;
    int degreey = 8;

    GlobalSliceMM global("../../../data/Bolshoi_7_clean_256.txt",
                         "/Users/auggiemarignier/Documents/PhD/TDT/massmapping/data/tutorial_prior.txt",
                         degreex,
                         degreey,
                         seed,
                         kmax,
                         wavelet_xy);
    BirthSlice birth(global);
    DeathSlice death(global);
    ValueSlice value(global);

    global.current_likelihood = global.likelihood(global.current_log_normalization);
    printf("Initial Likelihood: %f\n", global.current_likelihood);

    int *khistogram = new int[kmax];
    for (int i = 0; i < kmax; i++)
    {
        khistogram[i] = 0;
    }

    FILE *fp_ch = NULL;
    if (chain_history_initialise(global.ch,
                                 wavetree2d_sub_get_S_v(global.wt),
                                 global.current_likelihood,
                                 global.temperature,
                                 global.hierarchical->getparameter(0)) < 0)
    {
        fprintf(stderr, "error: failed to initialise chain history\n");
        return -1;
    }

    for (int i = 0; i < total; i++)
    {
        double u = global.random.uniform();
        std::cout << "i=" << i << "\t u=" << u << "\n";
    }
    return 0;
}