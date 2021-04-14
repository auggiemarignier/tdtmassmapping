#include <stdio.h>
#include <vector>

class mmobservations
{
public:
    double single_frequency_likelihood(size_t frequency_index,
                                       const double *model,
                                       double offset,
                                       double *residuals,
                                       double *residuals_normed,
                                       double &log_normalization);

    std::vector<double> single_frequency_predictions(std::vector<double> model);
};