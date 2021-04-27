#include <stdio.h>
#include <vector>

#include "hierarchicalmodel.hpp"

class mmobservations
{
public:
    double single_frequency_likelihood(
        std::vector<double> model,
        const hierarchicalmodel *hmodel,
        double *residuals,
        double *residuals_normed,
        double &log_normalization);

    std::vector<double> single_frequency_predictions(std::vector<double> model);
};