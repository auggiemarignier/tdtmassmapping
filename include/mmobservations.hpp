class mmobservations
{
public:
    double single_frequency_likelihood(size_t frequency_index,
                                       const double *model,
                                       double offset,
                                       const hierarchicalmodel *hmodel,
                                       double *residuals,
                                       double *residuals_normed,
                                       double &log_normalization);

    bool single_frequency_predictions(const double *model);
}