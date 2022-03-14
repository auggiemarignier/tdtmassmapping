#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <getopt.h>

#include "mmobservations.hpp"
#include "logging.hpp"
#include "utils.hpp"

static char short_options[] = "i:g:I:M:o:m:x:y:t:S:s:k:B:w:v:l:h:n:";
static struct option long_options[] = {
    {"input", required_argument, 0, 'i'},
    {"input_gamma", required_argument, 0, 'g'},
    {"output", required_argument, 0, 'o'},
    {"mask", required_argument, 0, 'm'},

    {"degree-x", required_argument, 0, 'x'},
    {"degree-y", required_argument, 0, 'y'},

    {"seed", required_argument, 0, 'S'},

    {"verbosity", required_argument, 0, 'v'},
    {"logfile", required_argument, 0, 'l'},
    {"help", no_argument, 0, 'h'},
    {"ngals", required_argument, 0, 'n'},

    {0, 0, 0, 0}};

static void usage(const char *pname);

int main(int argc, char *argv[])
{
    // Defaults
    char *input_kappa = nullptr;
    char *input_gamma = nullptr;
    char *output_prefix = nullptr;
    char *logfile = nullptr;
    char *maskfile = nullptr;

    int seed = 1;
    int verbosity = 1000;

    int degreex = 8;
    int degreey = 8;
    int super = 1;

    double *model;

    double ngal = 480.;     // galaxies per arcmin^2

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
            input_kappa = optarg;
            break;
        case 'g':
            input_gamma = optarg;
            break;
        case 'o':
            output_prefix = optarg;
            break;
        case 'm':
            maskfile = optarg;
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
        case 'S':
            seed = atoi(optarg);
            break;
        case 'v':
            verbosity = atoi(optarg);
            break;
        case 'l':
            logfile = optarg;
            break;
        case 'n':
            ngal = atof(optarg);
            break;
        case 'h':
        default:
            usage(argv[0]);
            return -1;
            break;
        }
    }
    Logger::open_log(logfile);

    // Check files
    if (input_kappa == NULL & input_gamma == NULL)
        throw ERROR("Please provide an input file\n");
    if (input_kappa != NULL & input_gamma != NULL)
        throw ERROR("Please provide only one of kappa or gamma file\n");
    if (output_prefix == NULL)
        throw ERROR("Please provide an output directory\n");

    mmobservations observations(1 << degreex, 1 << degreey, super);
    complexvector gamma;
    complexvector gamma_noisy;
    std::vector<double> covariance;

    if (input_kappa != NULL)
    {
        INFO("Reading in kappa map %s", input_kappa);
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
            throw ERROR("File not opened %s", input_kappa);
        }
        std::vector<int> mask(kappa.size(), 1);
        if (maskfile != nullptr)
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
        {
            mmobservations dummy_obs(1 << degreex, 1 << degreey, 1);
            gamma = dummy_obs.single_frequency_predictions(kappa);
        }
        const double sidelength = 10.; // arcmin
        INFO("%4.0f galaxies per pixel", ngal * std::pow(sidelength, 2) / kappa.size());
        bool aniso = false;
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
    else if (input_gamma != NULL)
    {
        INFO("Reading in gamma data %s", input_gamma);
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
            throw ERROR("File not opened %s", input_gamma);
        }
    }

    observations.set_observed_data(gamma_noisy);
    observations.set_sigmas(covariance);

    complexvector kappa_out(gamma_noisy.size());

    observations.kaiser_squires_inv(kappa_out, gamma_noisy);

    std::string filename = (input_kappa != NULL) ? "KS_kappa_synth.txt" : "KS_kappa_a520.txt";
    std::string output_file = mkfilename(output_prefix, filename.c_str());
    FILE *fp_out = fopen(output_file.c_str(), "w");
    if (fp_out == NULL)
    {
        fprintf(stderr, "error: failed to open mean file\n");
        return -1;
    }

    for (int j = 0; j < 1 << degreex; j++)
    {
        for (int i = 0; i < 1 << degreey; i++)
        {
            fprintf(fp_out, "%10.6f ", kappa_out[j * (1 << degreey) + i].real());
        }
        fprintf(fp_out, "\n");
    }

    fclose(fp_out);

    INFO("DONE");
}

static void usage(const char *pname)
{
    fprintf(stderr,
            "usage: %s [options]\n"
            "where options is one or more of:\n"
            "\n"
            " -i|--input <file>               Input kappa file.  One column\n"
            " -g|--input_gamma <file>         Input gamma file. Three columns (g1, g2, ngal)\n"
            " -m|--mask <file>                Mask file\n"
            " -o|--output <path>              Output prefix for output files\n"
            "\n"
            " -x|--degree-x <int>             Number of samples in x direction as power of 2\n"
            " -y|--degree-y <int>             Number of samples in y direction as power of 2\n"
            "\n"
            " -S|--seed <int>                 Random number seed\n"
            "\n"
            " -v|--verbosity <int>            Number steps between status printouts (0 = disable)\n"
            " -h|--help                       Show usage information\n"
            "\n",
            pname);
}