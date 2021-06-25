#include <iostream>
#include <gsl/gsl_rng.h>
#include <getopt.h>

#include "proposals.hpp"
#include "mmobservations.hpp"
#include "utils.hpp"

extern "C"
{
#include "slog.h"
};

static char short_options[] = "i:M:o:x:y:t:S:k:B:w:vh";
static struct option long_options[] = {
    {"input", required_argument, 0, 'i'},
    {"prior-file", required_argument, 0, 'M'},
    {"output", required_argument, 0, 'o'},

    {"degree-x", required_argument, 0, 'x'},
    {"degree-y", required_argument, 0, 'y'},

    {"total", required_argument, 0, 't'},
    {"seed", required_argument, 0, 'S'},

    {"kmax", required_argument, 0, 'k'},
    {"birth-probability", required_argument, 0, 'B'},

    {"wavelet-lateral", required_argument, 0, 'w'},

    {"verbosity", required_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},

    {0, 0, 0, 0}};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
    // Defaults
    char *input_obs = nullptr;
    char *prior_file = nullptr;
    char *output_prefix = nullptr;

    int total = 10000;
    int seed = 1;
    int verbosity = 1000;

    int kmax = 100;

    double Pb = 0.05;

    int wavelet_xy = 4;
    int degreex = 8;
    int degreey = 8;

    // Command line options
    int option_index = 0;
    while (true)
    {
        int c = getopt_long(argc, argv, short_options, long_options, &option_index);
        if (c == -1)
        {
            break;
        }
        switch (c)
        {
        case 'i':
            input_obs = optarg;
            break;
        case 'M':
            prior_file = optarg;
            break;
        case 'o':
            output_prefix = optarg;
            break;
        case 'x':
            degreex = atoi(optarg);
            if (degreex < 1 || degreex > 16)
            {
                fprintf(stderr, "error: degree x must be between 1 and 16 inclusive\n");
                return -1;
            }
            break;
        case 'y':
            degreey = atoi(optarg);
            if (degreey < 1 || degreey > 16)
            {
                fprintf(stderr, "error: degree y must be between 1 and 16 inclusive\n");
                return -1;
            }
            break;
        case 't':
            total = atoi(optarg);
            if (total < 1)
            {
                fprintf(stderr, "error: total must be greater than 0\n");
                return -1;
            }
            break;
        case 'S':
            seed = atoi(optarg);
            break;
        case 'k':
            kmax = atoi(optarg);
            if (kmax < 1)
            {
                fprintf(stderr, "error: kmax must be greater than 0\n");
                return -1;
            }
            break;
        case 'B':
            Pb = atof(optarg);
            if (Pb < 0.0 || Pb > 0.5)
            {
                fprintf(stderr, "error: birth probability must be between 0 and 0.5\n");
                return -1;
            }
            break;
        case 'w':
            wavelet_xy = atoi(optarg);
            if (wavelet_xy < 0 || wavelet_xy > GlobalProposal::WAVELET_MAX)
            {
                fprintf(stderr, "error: horizontal wavelet must be in range 0 .. %d\n", (int)GlobalProposal::WAVELET_MAX);
                return -1;
            }
            break;
        case 'v':
            verbosity = atoi(optarg);
            break;
        case 'h':
        default:
            usage(argv[0]);
            return -1;
            break;
        }
    }

    // Check files
    if (input_obs == NULL)
    {
        fprintf(stderr, "Please provide an input file\n");
        return -1;
    }
    if (prior_file == NULL)
    {
        fprintf(stderr, "Please provide a prior file\n");
        return -1;
    }
    if (output_prefix == NULL)
    {
        fprintf(stderr, "Please provide an output directory\n");
        return -1;
    }

    // Setup
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

    INFO("MCMC loop complete. Tidying up.");

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
    fprintf(fp, "%s", birth.write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, "%s", death.write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, "%s", value.write_long_stats().c_str());
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

    INFO("DONE");
    return 0;
}

static void usage(const char *pname)
{
    fprintf(stderr,
            "usage: %s [options]\n"
            "where options is one or more of:\n"
            "\n"
            " -i|--input <file>               Input observations file\n"
            " -M|--prior <file>               Prior/Proposal file\n"
            " -o|--output <path>              Output prefix for output files\n"
            "\n"
            " -x|--degree-x <int>             Number of samples in x direction as power of 2\n"
            " -y|--degree-y <int>             Number of samples in y direction as power of 2\n"
            "\n"
            " -t|--total <int>                Total number of iterations\n"
            " -S|--seed <int>                 Random number seed\n"
            "\n"
            " -k|--kmax <int>                 Max. no. of coefficients\n"
            "\n"
            " -B|--birth-probability <float>  Birth probability\n"
            "\n"
            " -w|--wavelet-xy <int>           Wavelet basis to use for x/y plane\n"
            "\n"
            " -v|--verbosity <int>            Number steps between status printouts (0 = disable)\n"
            " -h|--help                       Show usage information\n"
            "\n",
            pname);
}