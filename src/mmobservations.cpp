#include <vector>

#include "mmobservations.hpp"
extern "C"
{
#include "slog.h"
};

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
}