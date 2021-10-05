#include "utils.hpp"
#include <random>
#include <time.h>

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

std::complex<double> vector_mean(std::vector<std::complex<double>> &vec)
{
    std::complex<double> mean(0, 0);
    for (auto v : vec)
        mean += v;
    mean /= vec.size();
    return mean;
}

double vector_stddev(std::vector<std::complex<double>> &vec)
{
    double var = 0;
    std::complex<double> mean = vector_mean(vec);
    for (auto v : vec)
    {
        double absdiff = std::abs(v - mean);
        var += absdiff * absdiff;
    }
    var /= vec.size();
    return sqrt(var);
}

std::tuple<complexvector, std::vector<double>> add_gaussian_noise(const complexvector &input, const int &ngal, const int &sidelength)
{
    const int N = input.size();
    const double sigma_e = 0.37;               // intrinsic ellipticity dispersion
    const auto gals_pix = ngal * std::pow(sidelength, 2) / static_cast<double>(N);
    const auto A = sigma_e / std::sqrt(2.0 * gals_pix);

    // Initialize seed for Gaussian random number
    auto const seed = time(0);
    std::srand((unsigned int)seed);
    std::mt19937 mersenne(time(0));
    std::normal_distribution<> gaussian_dist(0, 1.0);

    // Create empty maps for output + covariance and initialize factors.
    complexvector output(N);
    std::vector<double> covariance(N);
    for (int i = 0; i < N; i++)
    {
        covariance[i] = A;
        output[i] = input[i] + A * std::complex<double>(gaussian_dist(mersenne), gaussian_dist(mersenne));
    }

    return std::make_tuple(output, covariance);
}
