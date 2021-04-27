#include <vector>
#include <iostream>

#include "mmobservations.hpp"
extern "C"
{
#include "slog.h"
};

mmobservations::mmobservations(
    std::vector<double> _obs,
    std::vector<double> _sigma) : obs(_obs),
                                  sigma(_sigma),
                                  n_obs(_obs.size()) {}

mmobservations::mmobservations(
    std::vector<double> _obs,
    double _sigma) : obs(_obs),
                     n_obs(_obs.size())
{
    for (int i = 0; i < n_obs; i++)
    {
        sigma[i] = _sigma;
    }
}

std::vector<double> mmobservations::single_frequency_predictions(
    std::vector<double> model)
{
    std::vector<double> predictions = model;
    return predictions;
}

double mmobservations::single_frequency_likelihood(
    std::vector<double> model,
    const hierarchicalmodel *hmodel,
    double *residuals,
    double *residuals_normed,
    double &log_normalization)
{
    std::vector<double> predictions = single_frequency_predictions(model);

    for (int i = 0; i < predictions.size(); i++)
    {
        residuals[i] = model[i] - predictions[i];
    }
    double loglikelihood = hmodel->nll(
        residuals,
        &sigma[0],
        n_obs,
        residuals_normed,
        log_normalization);

    return loglikelihood;
}