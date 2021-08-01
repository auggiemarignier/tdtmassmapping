#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "mmobservations.hpp"
#include "logging.hpp"

Observations::Observations(
    std::vector<double> _obs,
    std::vector<double> _sigma)
    : obs(_obs),
      sigma(_sigma),
      n_obs(_obs.size())
{
}

Observations::Observations(
    std::vector<double> _obs,
    double _sigma)
    : obs(_obs),
      n_obs(_obs.size())
{
    for (size_t i = 0; i < n_obs; i++)
    {
        sigma.push_back(_sigma);
    }
}

Observations::Observations(const char *filename)
{
    INFO("Opening file %s \n", filename);

    std::ifstream file(filename);
    double element;
    double mean = 0;
    double var = 0;
    double stddev = 0;
    if (file.is_open())
    {
        while (file >> element)
        {
            obs.push_back(element);
            mean += element;
        }
        n_obs = obs.size();
        mean /= n_obs;
        for (size_t i = 0; i < n_obs; i++)
        {
            var += (obs[i] - mean) * (obs[i] - mean);
        }
        var /= n_obs;
        stddev = sqrt(var);
        for (size_t i = 0; i < n_obs; i++)
        {
            sigma.push_back(stddev);
        }
        file.close();
    }
    else
    {
        throw ERROR("File not opened %s", filename);
    }
}

bool Observations::save_residuals(const char *filename,
                                  const double *residuals,
                                  const double *residuals_normed)
{
    const double *res = residuals;
    const double *resn = residuals_normed;
    FILE *fp = fopen(filename, "w");
    if (fp == NULL)
    {
        return false;
    }

    for (size_t i = 0; i < n_obs; i++)
    {
        fprintf(fp, "%15.9f %15.9f\n", *res, *resn);
        res++;
        resn++;
    }

    fclose(fp);
    return true;
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

    for (size_t i = 0; i < obs.size(); i++)
    {
        residuals[i] = obs[i] - predictions[i];
    }
    double loglikelihood = hmodel->nll(
        residuals,
        &sigma[0],
        n_obs,
        residuals_normed,
        log_normalization);

    return loglikelihood;
}