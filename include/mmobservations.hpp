#include <stdio.h>
#include <vector>

#include "hierarchicalmodel.hpp"

class mmobservations
{
public:
    // Constructor that reads a file
    // TODO: implement this
    mmobservations(const char *filename);

    // Constructor that takes in vectors
    mmobservations(std::vector<double> obs, std::vector<double> sigma);

    double single_frequency_likelihood(
        std::vector<double> model,
        const hierarchicalmodel *hmodel,
        double *residuals,
        double *residuals_normed,
        double &log_normalization);

    std::vector<double> single_frequency_predictions(std::vector<double> model);
};