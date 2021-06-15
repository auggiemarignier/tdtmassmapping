#include <iostream>
#include <gsl/gsl_rng.h>

#include "wavetomo2dutil.hpp"

#include "proposals.hpp"
#include "mmobservations.hpp"

extern "C"
{
#include "slog.h"
};

int main()
{
    const char *input_obs = "Bolshoi_7_clean_256.txt";
    const char *prior_file = "tutorial_prior.txt";
    const char *output_prefix = "../outputs/";

    int total = 10000;
    int seed = 1;
    int verbosity = 1000;

    int kmax = 100;

    double Pb = 0.05;

    int wavelet_xy = 4;
    int degreex = 8;
    int degreey = 8;

    mmobservations observations(input_obs);

    GlobalProposal global(&observations,
                          NULL,
                          prior_file,
                          degreex,
                          degreey,
                          seed,
                          kmax,
                          wavelet_xy);
    BirthProposal birth(global);
    DeathProposal death(global);
    ValueProposal value(global);

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

    for (int i = 0; i < total; i++) // start MCMC loop
    {
        double u = global.random.uniform();

        if (u < Pb) // Birth
        {
            if (birth.step() < 0)
            {
                fprintf(stderr, "error: failed to do birth step\n");
                return -1;
            }
        }
        else if (u < (2.0 * Pb)) // Death
        {
            if (death.step() < 0)
            {
                fprintf(stderr, "error: failed to do death step\n");
                return -1;
            }
        }
        else // Value
        {
            if (value.step() < 0)
            {
                fprintf(stderr, "error: failed to do value step\n");
                return -1;
            }
        }

        if (chain_history_full(global.ch))
        {
            // Flush chain history to file
            if (chain_history_write(global.ch,
                                    (ch_write_t)fwrite,
                                    fp_ch) < 0)
            {
                fprintf(stderr, "error: failed to write chain history segment to file\n");
                return -1;
            }
            if (chain_history_reset(global.ch) < 0)
            {
                fprintf(stderr, "error: failed to reset chain history\n");
                return -1;
            }
        }

        chain_history_change_t step;
        if (wavetree2d_sub_get_last_perturbation(global.wt, &step) < 0)
        {
            fprintf(stderr, "error: failed to get last step\n");
            return -1;
        }
        step.header.likelihood = global.current_likelihood;
        step.header.temperature = global.temperature;
        step.header.hierarchical = global.hierarchical->getparameter(0);
        if (chain_history_add_step(global.ch, &step) < 0)
        {
            fprintf(stderr, "error: failed to add step to chain history\n");
            return -1;
        }

        int current_k = wavetree2d_sub_coeff_count(global.wt);
        if (verbosity > 0 && (i + 1) % verbosity == 0)
        {
            INFO("%6d: %f %d dc %f lambda %f:\n",
                 i + 1,
                 global.current_likelihood,
                 current_k,
                 wavetree2d_sub_dc(global.wt),
                 global.hierarchical->getparameter(0));
            INFO(birth.write_long_stats().c_str());
            INFO(death.write_long_stats().c_str());
            INFO(value.write_long_stats().c_str());
        }
        khistogram[current_k - 1]++;
    } // end MCMC loop

    filename = mkfilename(output_prefix, "khistogram.txt");
    FILE *fp = fopen(filename.c_str(), "w");
    if (fp == NULL)
    {
        fprintf(stderr, "error: failed to create khistogram file\n");
        return -1;
    }
    for (int i = 0; i < kmax; i++)
    {
        fprintf(fp, "%d %d\n", i + 1, khistogram[i]);
    }
    fclose(fp);
    delete[] khistogram;

    // If there are remaining steps to save
    if (chain_history_nsteps(global.ch) > 1)
    {
        // Flush chain history to file
        if (chain_history_write(global.ch,
                                (ch_write_t)fwrite,
                                fp_ch) < 0)
        {
            fprintf(stderr, "error: failed to write chain history segment to file\n");
            return -1;
        }
    }
    fclose(fp_ch);
    chain_history_destroy(global.ch);

    filename = mkfilename(output_prefix, "acceptance.txt");
    fp = fopen(filename.c_str(), "w");
    if (fp == NULL)
    {
        fprintf(stderr, "error: failed to create acceptance file\n");
        return -1;
    }
    fprintf(fp, birth.write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, death.write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, value.write_long_stats().c_str());
    fprintf(fp, "\n");
    fclose(fp);

    filename = mkfilename(output_prefix, "final_model.txt");
    if (wavetree2d_sub_save(global.wt, filename.c_str()) < 0)
    {
        fprintf(stderr, "error: failed to save final model\n");
        return -1;
    }

    filename = mkfilename(output_prefix, "residuals.txt");
    if (!global.save_residuals(filename.c_str()))
    {
        ERROR("Failed to create residuals file");
        return -1;
    }

    filename = mkfilename(output_prefix, "residuals_hist.txt");
    if (!global.save_residual_histogram(filename.c_str()))
    {
        ERROR("Failed to save residual histogram");
        return -1;
    }

    filename = mkfilename(output_prefix, "residuals_cov.txt");
    if (!global.save_residual_covariance(filename.c_str()))
    {
        ERROR("Failed to save residual covariance");
        return -1;
    }

    return 0;
}