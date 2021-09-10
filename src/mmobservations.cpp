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

std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>> mmobservations::init_fft_2d(const uint &imsizey, const uint &imsizex)
{
    uint imsize = imsizex * imsizey;
    fftw_complex *in, *out;

    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    auto del = [](fftw_plan_s *plan)
    { fftw_destroy_plan(plan); };
    std::shared_ptr<fftw_plan_s> plan_forward(fftw_plan_dft_2d(imsizey, imsizex, in, out, FFTW_FORWARD, FFTW_MEASURE), del);
    std::shared_ptr<fftw_plan_s> plan_inverse(fftw_plan_dft_2d(imsizey, imsizex, in, out, FFTW_BACKWARD, FFTW_MEASURE), del);

    auto forward = [=](fftw_complex *output, const fftw_complex *input)
    {
        fftw_execute_dft(plan_forward.get(), const_cast<fftw_complex *>(input), output);
        for (int i = 0; i < (int)(imsize); i++)
        {
            output[i][0] /= std::sqrt(imsize);
            output[i][1] /= std::sqrt(imsize);
        }
    };

    auto backward = [=](fftw_complex *output, const fftw_complex *input)
    {
        fftw_execute_dft(plan_inverse.get(), const_cast<fftw_complex *>(input), output);
        for (int i = 0; i < (int)(imsize); i++)
        {
            output[i][0] /= std::sqrt(imsize);
            output[i][1] /= std::sqrt(imsize);
        }
    };

    return std::make_tuple(forward, backward);
}

std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>> mmobservations::build_lensing_kernels(const uint &imsizey, const uint &imsizex)
{
    const uint imsize = imsizex * imsizey;
    const int n = (int)std::sqrt(imsize);
    double kx, ky;

    complexvector lensing_kernel;
    complexvector adjoint_kernel;
    lensing_kernel.reserve(imsize);
    adjoint_kernel.reserve(imsize);

    for (int i = 0; i < n; i++)
    {
        ky = (2 * i - 2 < n) ? (static_cast<double>(i)) : (static_cast<double>(i - n));

        for (int j = 0; j < n; j++)
        {
            kx = (2 * j - 2 < n) ? (static_cast<double>(j)) : (static_cast<double>(j - n));

            if ((kx != 0.0) || (ky != 0.0))
            {
                double real = (ky * ky - kx * kx) / (kx * kx + ky * ky);
                double imag = (2.0 * kx * ky) / (kx * kx + ky * ky);
                lensing_kernel.emplace_back(real, imag);
                adjoint_kernel.emplace_back(real, -imag);
            }
            else
            {
                lensing_kernel.emplace_back(0, 0);
                adjoint_kernel.emplace_back(0, 0);
            }
        }
    }

    auto forward = [=](fftw_complex *output, const fftw_complex *input)
    {
        for (int i = 0; i < (int)imsize; i++)
        {
            output[i][0] = lensing_kernel[i].real() * input[i][0];
            output[i][1] = lensing_kernel[i].imag() * input[i][1];
        }
    };

    auto adjoint = [=](fftw_complex *output, const fftw_complex *input)
    {
        for (int i = 0; i < (int)imsize; i++)
        {
            output[i][0] = adjoint_kernel[i].real() * input[i][0];
            output[i][1] = adjoint_kernel[i].imag() * input[i][1];
        }
    };

    return std::make_tuple(forward, adjoint);
}