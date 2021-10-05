#pragma once

#include "wavetree.h"
#include <string>
#include <complex>
#include <vector>
#include <math.h>

typedef std::vector<std::complex<double>> complexvector;

std::string enum_to_string(wavetree_perturb_t type);

std::string mkfilename(const char *prefix, const char *file);

std::string mkformatstring(const char *fmt, ...);

std::complex<double> vector_mean(complexvector &v);
double vector_stddev(complexvector &v);

std::tuple<complexvector, std::vector<double>> add_gaussian_noise(const complexvector &input, const int &ngal, const int sidelegnth);

namespace statistics {
std::tuple<double, double> run_statistics(const std::vector<double>& truth, const std::vector<double>& estimate);
double snr(const std::vector<double>& truth, const std::vector<double>& estimate); 
double pearson_correlation(const std::vector<double>& truth, const std::vector<double>& estimate);
};  // namespace statistics