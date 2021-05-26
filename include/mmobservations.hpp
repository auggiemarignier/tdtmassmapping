#pragma once

#include <stdio.h>
#include <vector>

#include "hierarchicalmodel.hpp"

class mmobservations
{
public:
    // Constructor that takes in vectors
    mmobservations(std::vector<double> _obs, std::vector<double> _sigma);

    // Constructor that takes a vector of obs and the stddev
    mmobservations(std::vector<double> _obs, double _sigma);

    double single_frequency_likelihood(
        std::vector<double> model,
        const hierarchicalmodel *hmodel,
        double *residuals,
        double *residuals_normed,
        double &log_normalization);

    std::vector<double> single_frequency_predictions(std::vector<double> model);

    std::vector<double> obs;
    std::vector<double> sigma;
    size_t n_obs;
};