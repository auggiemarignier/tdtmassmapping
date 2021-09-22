#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>

#include "mmobservations.hpp"
#include "logging.hpp"

#if 0
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
#endif

mmobservations::mmobservations(const uint _imsizex, const uint _imsizey)
    : Observations(),
      imsizex(_imsizex),
      imsizey(_imsizey),
      imsize(_imsizey * _imsizex)
{
    auto fft_tuple = init_fft_2d();
    fft = std::get<0>(fft_tuple);
    ifft = std::get<1>(fft_tuple);

    auto operator_tuple = build_lensing_kernels();
    D = std::get<0>(operator_tuple);
    Dinv = std::get<1>(operator_tuple);
    Dadj = std::get<2>(operator_tuple);
};

complexvector mmobservations::single_frequency_predictions(
    complexvector model)
{
    complexvector predictions = model;
    return predictions;
}

double mmobservations::single_frequency_likelihood(
    complexvector model,
    std::complex<double> *residuals,
    std::complex<double> *residuals_normed,
    double &log_normalization)
{
    complexvector predictions = single_frequency_predictions(model);

    double loglikelihood = 0.0;
    for (size_t i = 0; i < imsize; i++)
    {
        residuals[i] = obs[i] - predictions[i];
        residuals_normed[i] = residuals[i] / sigma[i];
        loglikelihood += std::abs(residuals_normed[i] * residuals_normed[i]) * 0.5;
        log_normalization += log(sigma[i]);
    }

    return loglikelihood;
}

std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>> mmobservations::init_fft_2d()
{
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

std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>> mmobservations::build_lensing_kernels()
{
    const int n = (int)std::sqrt(imsize);
    double kx, ky;

    complexvector lensing_kernel;
    complexvector adjoint_kernel;
    lensing_kernel.reserve(imsize);
    adjoint_kernel.reserve(imsize);

    for (int i = 0; i < n; i++)
    {
        ky = (2 * i < n) ? (static_cast<double>(i)) : (static_cast<double>(i - n));

        for (int j = 0; j < n; j++)
        {
            kx = (2 * j < n) ? (static_cast<double>(j)) : (static_cast<double>(j - n));

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
            double lkr = lensing_kernel[i].real();
            double lki = lensing_kernel[i].imag();
            double inr = input[i][0];
            double ini = input[i][1];

            output[i][0] = lkr * inr - lki * ini;
            output[i][1] = lkr * ini + lki * inr;
        }
    };

    auto inverse = [=](fftw_complex *output, const fftw_complex *input)
    {
        for (int i = 0; i < (int)imsize; i++)
        {
            double lkr = lensing_kernel[i].real();
            double lki = lensing_kernel[i].imag();
            double norm = lkr * lkr + lki * lki;
            double inr = input[i][0];
            double ini = input[i][1];

            output[i][0] = (i == 0) ? 0 : (inr * lkr + ini * lki) / norm;
            output[i][1] = (i == 0) ? 0 : (ini * lkr - inr * lki) / norm;
        }
    };

    auto adjoint = [=](fftw_complex *output, const fftw_complex *input)
    {
        for (int i = 0; i < (int)imsize; i++)
        {
            double lkr = adjoint_kernel[i].real();
            double lki = adjoint_kernel[i].imag();
            double inr = input[i][0];
            double ini = input[i][1];

            output[i][0] = lkr * inr - lki * ini;
            output[i][1] = lkr * ini + lki * inr;
        }
    };

    return std::make_tuple(forward, inverse, adjoint);
}

void mmobservations::kaiser_squires(fftw_complex *output, const fftw_complex *input)
{
    fftw_complex *temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *temp2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fft(temp, input);
    D(temp2, temp);
    ifft(output, temp2);
}

void mmobservations::kaiser_squires_inv(fftw_complex *output, const fftw_complex *input)
{
    fftw_complex *temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *temp2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fft(temp, input);
    Dinv(temp2, temp);
    ifft(output, temp2);
}

void mmobservations::kaiser_squires_adj(fftw_complex *output, const fftw_complex *input)
{
    fftw_complex *temp = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *temp2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fft(temp, input);
    Dadj(temp2, temp);
    ifft(output, temp2);
}

void mmobservations::set_observed_data(complexvector &_obs)
{
    if (_obs.size() == imsize)
        obs = _obs;
    else
        ERROR("Input data has size %i.  Expected size %i.", _obs.size(), imsize);
}

void mmobservations::set_sigmas(std::vector<double> &_sigmas)
{
    if (_sigmas.size() == imsize)
        sigma = _sigmas;
    else
        ERROR("Input data has size %i.  Expected size %i.", _sigmas.size(), imsize);
}