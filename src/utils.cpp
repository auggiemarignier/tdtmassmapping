#include "utils.hpp"

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
    for (auto v:vec)
    {
        double absdiff = std::abs(v - mean);
        var += absdiff * absdiff;
    }
    var /= vec.size();
    return sqrt(var);
}
