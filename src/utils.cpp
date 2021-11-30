#include "utils.hpp"
#include <random>
#include <time.h>
#include "logging.hpp"

std::string enum_to_string(wavetree_perturb_t type)
{
    switch (type)
    {
    case WT_PERTURB_INVALID:
        return "Invalid";
    case WT_PERTURB_NONE:
        return "None";
    case WT_PERTURB_BIRTH:
        return "Birth";
    case WT_PERTURB_DEATH:
        return "Death";
    case WT_PERTURB_VALUE:
        return "Value";
    case WT_PERTURB_MOVE:
        return "Move";
    case WT_PERTURB_HIERARCHICAL:
        return "Hierarchical";
    case WT_PERTURB_PTEXCHANGE:
        return "PT Exchange";
    case WT_PERTURB_PTMODELEXCHANGE:
        return "PT Model Exchange";
    }
}

std::string mkfilename(const char *prefix, const char *file)
{
    if (prefix == nullptr)
    {
        return std::string(file);
    }
    else
    {
        return std::string(prefix) + file;
    }
}

std::string mkformatstring(const char *fmt, ...)
{
    static char *buffer = nullptr;
    static int buffer_size = -1;

    if (buffer == nullptr)
    {
        buffer_size = 512;
        buffer = new char[buffer_size];
    }

    va_list ap;
    int size;

    va_start(ap, fmt);
    size = vsnprintf(buffer, buffer_size, fmt, ap);
    while (size >= buffer_size)
    {
        delete[] buffer;
        buffer_size *= 2;
        buffer = new char[buffer_size];
        size = vsnprintf(buffer, buffer_size, fmt, ap);
    }
    va_end(ap);

    return std::string(buffer);
}

std::tuple<complexvector, std::vector<double>> add_gaussian_noise(const complexvector &input, const double &ngal, const double &sidelength, const bool aniso)
{
    const int N = input.size();
    const double sigma_e = 0.37; // intrinsic ellipticity dispersion
    std::vector<double> ngal_an;
    double gals_pix;
    double A;
    // Initialize seed for Gaussian random number
    std::mt19937 mersenne;
    std::normal_distribution<> gaussian_dist(0, 1.0);

    // Create empty maps for output + covariance and initialize factors.
    complexvector output(N);
    std::vector<double> covariance(N);

    for (int i = 0; i < N; i++)
    {
        if (aniso)
            ngal_an.push_back(ngal * std::abs(input[i]) + 1);
        else
            ngal_an.push_back(ngal);

        gals_pix = ngal_an[i] * std::pow(sidelength, 2) / static_cast<double>(N);
        A = sigma_e / std::sqrt(2.0 * gals_pix);
        covariance[i] = A;
        output[i] = input[i] + A * std::complex<double>(gaussian_dist(mersenne), gaussian_dist(mersenne));
    }

    return std::make_tuple(output, covariance);
}

namespace statistics
{
    std::tuple<double, double> run_statistics(const std::vector<double> &truth,
                                              const std::vector<double> &estimate)
    {
        return std::make_tuple(snr(truth, estimate), pearson_correlation(truth, estimate));
    };

    double snr(const std::vector<double> &truth, const std::vector<double> &estimate)
    {
        if (truth.size() != estimate.size())
        {
            WARNING("Vectors are of different sizes");
            return -1;
        }
        else
        {
            double l2_true = 0.;
            double l2_diff_to_true = 0.;
            for (int i = 0; i < truth.size(); i++)
            {
                l2_true += truth[i] * truth[i];
                l2_diff_to_true += std::pow((estimate[i] - truth[i]), 2);
            }
            return 10.0 * std::log10(l2_true / l2_diff_to_true);
        };
    }

    double snr(const complexvector &truth, const complexvector &estimate)
    {
        if (truth.size() != estimate.size())
        {
            WARNING("Vectors are of different sizes");
            return -1;
        }
        else
        {
            double l2_true = 0.;
            double l2_diff_to_true = 0.;
            for (int i = 0; i < truth.size(); i++)
            {
                l2_true += std::pow(std::abs(truth[i]), 2);
                l2_diff_to_true += std::pow(std::abs(estimate[i] - truth[i]), 2);
            }
            return 10.0 * std::log10(l2_true / l2_diff_to_true);
        };
    }

    double pearson_correlation(const std::vector<double> &truth, const std::vector<double> &estimate)
    {
        if (truth.size() != estimate.size())
        {
            WARNING("Vectors are of different sizes");
            return -1;
        }
        else
        {
            double numerator = 0.0;
            double sum_truth_squares = 0.0;
            double sum_estimate_squares = 0.0;
            const double mean_estimate = vector_mean(estimate);
            const double mean_truth = vector_mean(truth);
            for (uint i = 0; i < truth.size(); i++)
            {
                numerator += (truth[i] - mean_truth) * (estimate[i] - mean_estimate);
                sum_truth_squares += std::pow((truth[i] - mean_truth), 2);
                sum_estimate_squares += std::pow((estimate[i] - mean_estimate), 2);
            }
            return numerator / (std::sqrt(sum_truth_squares) * std::sqrt(sum_estimate_squares));
        }
    };
} // namespace statistics