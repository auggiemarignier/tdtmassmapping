#pragma once

#include <stdio.h>
#include <vector>

#include "hierarchicalmodel.hpp"

class Observations
{
    virtual double single_frequency_likelihood(std::vector<double> model,
        const hierarchicalmodel *hmodel,
        double *residuals,
        double *residuals_normed,
        double &log_normalization) = 0;

    virtual std::vector<double> single_frequency_predictions(std::vector<double> model) = 0;
};

class mmobservations : Observations
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
        double &log_normalization) override;

    std::vector<double> single_frequency_predictions(std::vector<double> model) override;

    std::vector<double> obs;
    std::vector<double> sigma;
    size_t n_obs;
};