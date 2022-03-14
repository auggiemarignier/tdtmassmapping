#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <getopt.h>
#include <nlohmann/json.hpp>

#include "proposals.hpp"
#include "mmobservations.hpp"
#include "utils.hpp"
#include "globalprop.cpp"
#include "logging.hpp"

static void usage(const char *pname);

using json = nlohmann::json;

int main(int argc, char *argv[])
{
    // Defaults
    int total = 10000;
    int seed = 1;
    int verbosity = 1000;

    int kmax = 100;

    double Pb = 0.05;

    int wavelet_xy = 4;
    int degreex = 8;
    int degreey = 8;
    int super = 1;

    double ngal = 480.; // galaxies per arcmin^2

    double *model;

    // Load configuration
    json j;
    std::fstream cfile("data/config.json");
    cfile >> j;

    Logger::open_log(j["inputs"]["logfile"].get<json::string_t>());
    std::string input_kappa = j["WL_simulations"]["input_kappa"].get<json::string_t>();
    std::string maskfile = j["WL_simulations"]["maskfile"].get<json::string_t>();
    std::string input_gamma = j["inputs"]["input_gamma"].get<json::string_t>();
    const char *prior_file = const_cast<char *>(j["inputs"]["prior_file"].get<json::string_t>().c_str());
    const char *initial_model = const_cast<char *>(j["inputs"]["initial_model"].get<json::string_t>().c_str());
    const char *output_prefix = const_cast<char *>(j["inputs"]["output_prefix"].get<json::string_t>().c_str());
    LOG("%s", prior_file);

    // Check files
    if (input_kappa.empty() & input_gamma.empty())
        throw ERROR("Please provide an input file\n");
    if (!input_kappa.empty() & !input_gamma.empty())
        throw ERROR("Please provide only one of kappa or gamma file\n");
    if (prior_file == NULL)
        throw ERROR("Please provide a prior file\n");
    if (output_prefix == NULL)
        throw ERROR("Please provide an output directory\n");

    mmobservations observations(1 << degreex, 1 << degreey, super);
    complexvector gamma;
    complexvector gamma_noisy;
    std::vector<double> covariance;

    if (!input_kappa.empty())
    {
        INFO("Reading in kappa map %s", input_kappa.c_str());
        std::ifstream file(input_kappa);
        complexvector kappa;
        double element;
        if (file.is_open())
        {
            while (file >> element)
                kappa.emplace_back(element, 0);

            if (1 << degreex * 1 << degreey != kappa.size())
                throw ERROR("Incorrect image size");

            file.close();
        }
        else
        {
            throw ERROR("File not opened %s", input_kappa.c_str());
        }
        std::vector<int> mask(kappa.size(), 1);
        if (!maskfile.empty())
        {
            std::ifstream mfile(maskfile);
            int melement;
            int i = 0;
            if (mfile.is_open())
            {
                while (mfile >> melement)
                {
                    mask[i] = melement;
                    i++;
                }
            }
        }
        mmobservations dummy_obs(1 << degreex, 1 << degreey, 1);
        gamma = dummy_obs.single_frequency_predictions(kappa);
        const double sidelength = 8.; // arcmin
        INFO("%4.0f galaxies per pixel", ngal * std::pow(sidelength, 2) / kappa.size());
        bool aniso = true;
        auto noise_tuple = add_gaussian_noise(gamma, ngal, sidelength, aniso);
        gamma_noisy = std::get<0>(noise_tuple);
        covariance = std::get<1>(noise_tuple);
        for (int i = 0; i < gamma_noisy.size(); i++)
        {
            gamma_noisy[i] *= mask[i];
            covariance[i] /= mask[i]; // infinite covariance where mask=0
        }
        double data_snr = statistics::snr(gamma, gamma_noisy);
        INFO("Input data SNR %10.6f dB", data_snr);
    }
    else if (!input_gamma.empty())
    {
        INFO("Reading in gamma data %s", input_gamma.c_str());
        std::ifstream file(input_gamma);
        double re, im;
        int ngal;
        if (file.is_open())
        {
            while (file >> re >> im >> ngal)
            {
                gamma_noisy.emplace_back(re, im);
                covariance.emplace_back(0.37 * 1e-5 / sqrt(2. * ngal)); // infinite covariance where ngal=0
            }

            if (1 << degreex * 1 << degreey != gamma_noisy.size())
                throw ERROR("Incorrect image size");

            file.close();
        }
        else
        {
            throw ERROR("File not opened %s", input_gamma.c_str());
        }
    }

    observations.set_observed_data(gamma_noisy);
    observations.set_sigmas(covariance);

    GlobalProposal global(&observations,
                          initial_model,
                          prior_file,
                          degreex + super - 1,
                          degreey + super - 1,
                          seed,
                          kmax,
                          wavelet_xy);
    BirthProposal birth(global);
    DeathProposal death(global);
    ValueProposal value(global);

    global.current_likelihood = global.likelihood(global.current_log_normalization);
    INFO("Initial Likelihood: %f\n", global.current_likelihood);

    int *khistogram = new int[kmax];
    for (int i = 0; i < kmax; i++)
        khistogram[i] = 0;

    FILE *fp_ch = NULL;
    if (chain_history_initialise(global.ch,
                                 wavetree2d_sub_get_S_v(global.wt),
                                 global.current_likelihood,
                                 global.temperature,
                                 1.0) < 0)
        throw ERROR("failed to initialise chain history\n");

    std::string filename = mkfilename(output_prefix, "ch.dat");
    fp_ch = fopen(filename.c_str(), "w");
    if (fp_ch == NULL)
        throw ERROR("failed to create chain history file\n");

    Logger::flush();
    for (int i = 0; i < total; i++) // start MCMC loop
    {
        double u = global.random.uniform();

        if (u < Pb) // Birth
        {
            if (birth.step() < 0)
                throw ERROR("failed to do birth step\n");
        }
        else if (u < (2.0 * Pb)) // Death
        {
            if (death.step() < 0)
                throw ERROR("failed to do death step\n");
        }
        else // Value
        {
            if (value.step() < 0)
                throw ERROR("failed to do value step\n");
        }

        if (chain_history_full(global.ch))
        {
            // Flush chain history to file
            if (chain_history_write(global.ch,
                                    (ch_write_t)fwrite,
                                    fp_ch) < 0)
                throw ERROR("failed to write chain history segment to file\n");
            if (chain_history_reset(global.ch) < 0)
                throw ERROR("failed to reset chain history\n");
        }

        chain_history_change_t step;
        if (wavetree2d_sub_get_last_perturbation(global.wt, &step) < 0)
            throw ERROR("failed to get last step\n");

        step.header.likelihood = global.current_likelihood;
        step.header.temperature = global.temperature;
        step.header.hierarchical = 1.0;
        if (chain_history_add_step(global.ch, &step) < 0)
            throw ERROR("failed to add step to chain history\n");

        int current_k = wavetree2d_sub_coeff_count(global.wt);
        if (verbosity > 0 && (i + 1) % verbosity == 0)
        {
            INFO("Iteration %6d: Current Likelihood: %f Current k: %d\n",
                 i + 1,
                 global.current_likelihood,
                 current_k);
            INFO(birth.write_long_stats().c_str());
            INFO(death.write_long_stats().c_str());
            INFO(value.write_long_stats().c_str());
            LOG("\n");
            Logger::flush();
        }
        khistogram[current_k - 1]++;
    } // end MCMC loop

    INFO("MCMC loop complete. Tidying up.");

    filename = mkfilename(output_prefix, "khistogram.txt");
    FILE *fp = fopen(filename.c_str(), "w");
    if (fp == NULL)
        throw ERROR("failed to create khistogram file\n");

    for (int i = 0; i < kmax; i++)
        fprintf(fp, "%d %d\n", i + 1, khistogram[i]);

    fclose(fp);
    delete[] khistogram;

    // If there are remaining steps to save
    if (chain_history_nsteps(global.ch) > 1)
    {
        // Flush chain history to file
        if (chain_history_write(global.ch,
                                (ch_write_t)fwrite,
                                fp_ch) < 0)
            throw ERROR("failed to write chain history segment to file\n");
    }
    fclose(fp_ch);
    chain_history_destroy(global.ch);

    filename = mkfilename(output_prefix, "coeff_mean.txt");
    FILE *fp_mean = fopen(filename.c_str(), "w");
    filename = mkfilename(output_prefix, "coeff_std.txt");
    FILE *fp_std = fopen(filename.c_str(), "w");
    filename = mkfilename(output_prefix, "coeff_n.txt");
    FILE *fp_n = fopen(filename.c_str(), "w");
    int index = 0;
    double mean;
    double std;
    int n;
    for (int j = 0; j < global.height; j++)
    {
        for (int i = 0; i < global.width; i++)
        {
            n = coefficient_histogram_get_coefficient_mean_std(global.coeff_hist, index, &mean, &std);
            fprintf(fp_mean, "%10.6f ", mean);
            fprintf(fp_std, "%10.6f ", std);
            fprintf(fp_n, "%i ", n);
            index++;
        }
        fprintf(fp, "\n");
    }
    fclose(fp_mean);
    fclose(fp_std);
    fclose(fp_n);

    filename = mkfilename(output_prefix, "acceptance.txt");
    fp = fopen(filename.c_str(), "w");
    if (fp == NULL)
        throw ERROR("failed to create acceptance file\n");
    fprintf(fp, "%s", birth.write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, "%s", death.write_long_stats().c_str());
    fprintf(fp, "\n");
    fprintf(fp, "%s", value.write_long_stats().c_str());
    fprintf(fp, "\n");
    fclose(fp);

    filename = mkfilename(output_prefix, "final_model.txt");
    if (wavetree2d_sub_save(global.wt, filename.c_str()) < 0)
        throw ERROR("failed to save final model\n");

    filename = mkfilename(output_prefix, "final_model_pix.txt");
    model = new double[global.size];
    memset(model, 0, sizeof(double) * global.size);
    if (wavetree2d_sub_map_to_array(global.wt, model, global.size) < 0)
        throw ERROR("Failed to map model to array\n");

    if (generic_lift_inverse2d(model,
                               global.width,
                               global.height,
                               global.width,
                               global.workspace,
                               global.xywaveletf,
                               global.xywaveletf,
                               SUBTILE) < 0)
        throw ERROR("Failed to do inverse transform on coefficients\n");

    fp = fopen(filename.c_str(), "w");
    if (fp == NULL)
        throw ERROR("failed to create final_model_pix.txt file\n");

    for (int j = 0; j < global.height; j++)
    {
        for (int i = 0; i < global.width; i++)
        {
            fprintf(fp, "%10.6f ", model[j * global.width + i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    INFO("DONE");
    return 0;
}

static void usage(const char *pname)
{
    fprintf(stderr,
            "usage: %s [options]\n"
            "where options is one or more of:\n"
            "\n"
            " -i|--input <file>               Input kappa file.  One column\n"
            " -g|--input_gamma <file>         Input gamma file. Three columns (g1, g2, ngal)\n"
            " -I|--inital_model <file>        Initial model file for restarts\n"
            " -M|--prior <file>               Prior/Proposal file\n"
            " -m|--mask <file>                Mask file\n"
            " -o|--output <path>              Output prefix for output files\n"
            "\n"
            " -x|--degree-x <int>             Number of samples in x direction as power of 2\n"
            " -y|--degree-y <int>             Number of samples in y direction as power of 2\n"
            "\n"
            " -t|--total <int>                Total number of iterations\n"
            " -S|--seed <int>                 Random number seed\n"
            " -s|--super <int>                Super-resolution factor. Default = 1 i.e. no superresolution\n"
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