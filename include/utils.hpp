#pragma once

#include "wavetree.h"
#include <string>
#include <complex>
#include <vector>
#include <math.h>

std::string enum_to_string(wavetree_perturb_t type);

std::string mkfilename(const char *prefix, const char *file);

std::string mkformatstring(const char *fmt, ...);

std::complex<double> vector_mean(std::vector<std::complex<double>> &v);
double vector_stddev(std::vector<std::complex<double>> &v);