#include <iostream>
#include <gsl/gsl_rng.h>

#include "wavetomo2dutil.hpp"

#include "proposals.hpp"

int main()
{
    const char *input_obs = "Bolshoi_7_clean_256.txt";
    const char *prior_file = "tutorial_prior.txt";
    const char *output_prefix = "../outputs/";

    int total = 10;
    int seed = 1;

    int kmax = 100;

    double Pb = 0.05;

    int wavelet_xy = 4;
    int degreex = 8;
    int degreey = 8;

    GlobalSliceMM global(input_obs,
                         prior_file,
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

    std::string filename = mkfilename(output_prefix, "ch.dat");
    fp_ch = fopen(filename.c_str(), "w");
    if (fp_ch == NULL)
    {
        fprintf(stderr, "error: failed to create chain history file\n");
        return -1;
    }

    for (int i = 0; i < total; i++)
    {
        double u = global.random.uniform();

        if (u < Pb)
        {

            //
            // Birth
            //
            if (birth.step() < 0)
            {
                fprintf(stderr, "error: failed to do birth step\n");
                return -1;
            }
        }
        else if (u < (2.0 * Pb))
        {

            //
            // Death
            //
            if (death.step() < 0)
            {
                fprintf(stderr, "error: failed to do death step\n");
                return -1;
            }
        }
        else
        {

            //
            // Value
            //
            if (value.step() < 0)
            {
                fprintf(stderr, "error: failed to do value step\n");
                return -1;
            }
        }
    }

    fclose(fp_ch);

    return 0;
}