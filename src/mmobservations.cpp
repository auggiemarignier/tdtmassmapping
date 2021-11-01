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

double Observations::single_frequency_likelihood(
    complexvector &model,
    double &log_normalization)
{
    complexvector predictions = single_frequency_predictions(model);
    double loglikelihood = 0.0;
    for (size_t i = 0; i < n_obs; i++)
    {
        std::complex<double> res = obs[i] - predictions[i];
        std::complex<double> res_normed = res / sigma[i];
        loglikelihood += std::norm(res_normed) * 0.5;
        log_normalization += log(sigma[i]);
    }

    return loglikelihood;
}

void Observations::set_observed_data(complexvector &_obs)
{
    obs = _obs;
    n_obs = _obs.size();
}

void Observations::set_sigmas(std::vector<double> &_sigmas)
{
    if (_sigmas.size() == n_obs)
        sigma = _sigmas;
    else if (_sigmas.size() == 1)
    {
        sigma.reserve(n_obs);
        for (uint i = 0; i < n_obs; i++)
            sigma.emplace_back(_sigmas[0]);
    }
    else
        ERROR("Input data has size %i.  Expected size %i.", _sigmas.size(), n_obs);
}

complexvector Identity::single_frequency_predictions(complexvector &model)
{
    complexvector predictions = model;
    return predictions;
}

mmobservations::mmobservations(const uint _imsizex, const uint _imsizey, const uint _super)
    : Observations(),
      imsizex(_imsizex),
      imsizey(_imsizey),
      imsize(_imsizey * _imsizex),
      super(_super),
      superimsizex(_imsizex),
      superimsizey(_imsizey),
      superimsize(_imsizex * _imsizey)
{
    if (super > 1)
    {
        superimsizex = imsizex << (super - 1);
        superimsizey = imsizey << (super - 1);
        superimsize = superimsizex * superimsizey;
    }
    auto fft_tuple = init_fft_2d(imsizex, imsizey);
    fft = std::get<0>(fft_tuple);
    ifft = std::get<1>(fft_tuple);
    if (super > 1)
    {
        auto s_fft_tuple = init_fft_2d(superimsizex, superimsizey);
        s_fft = std::get<0>(s_fft_tuple);
        s_ifft = std::get<1>(s_fft_tuple);
    }

    auto operator_tuple = build_lensing_kernels();
    D = std::get<0>(operator_tuple);
    Dinv = std::get<1>(operator_tuple);
    Dadj = std::get<2>(operator_tuple);
};

complexvector mmobservations::single_frequency_predictions(complexvector &kappa)
{
    complexvector gamma(imsize);
    kaiser_squires(gamma, kappa);
    return gamma;
}

std::tuple<std::function<void(complexvector &, const complexvector &)>, std::function<void(complexvector &, const complexvector &)>> mmobservations::init_fft_2d(const uint _imsizex, const uint _imsizey)
{
    const uint _imsize = _imsizex * _imsizey;
    complexvector in(_imsize);
    complexvector out(_imsize);

    auto del = [](fftw_plan_s *plan)
    { fftw_destroy_plan(plan); };
    std::shared_ptr<fftw_plan_s> plan_forward(
        fftw_plan_dft_2d(
            _imsizey,
            _imsizex,
            reinterpret_cast<fftw_complex *>(&in[0]),
            reinterpret_cast<fftw_complex *>(&out[0]),
            FFTW_FORWARD,
            FFTW_MEASURE),
        del);
    std::shared_ptr<fftw_plan_s> plan_inverse(
        fftw_plan_dft_2d(
            _imsizey,
            _imsizex,
            reinterpret_cast<fftw_complex *>(&in[0]),
            reinterpret_cast<fftw_complex *>(&out[0]),
            FFTW_BACKWARD,
            FFTW_MEASURE),
        del);

    auto forward = [=](complexvector &output, const complexvector &input)
    {
        fftw_execute_dft(
            plan_forward.get(),
            const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(&input[0])),
            reinterpret_cast<fftw_complex *>(&output[0]));
        for (int i = 0; i < (int)(_imsize); i++)
            output[i] /= std::sqrt(_imsize);
    };

    auto backward = [=](complexvector &output, const complexvector &input)
    {
        fftw_execute_dft(
            plan_inverse.get(),
            const_cast<fftw_complex *>(reinterpret_cast<const fftw_complex *>(&input[0])),
            reinterpret_cast<fftw_complex *>(&output[0]));
        for (int i = 0; i < (int)(_imsize); i++)
            output[i] /= std::sqrt(_imsize);
    };

    return std::make_tuple(forward, backward);
}

std::tuple<std::function<void(complexvector &, const complexvector &)>, std::function<void(complexvector &, const complexvector &)>, std::function<void(complexvector &, const complexvector &)>> mmobservations::build_lensing_kernels()
{
    const int n = (int)std::sqrt(imsize);
    double kx, ky;

    complexvector lensing_kernel(imsize);
    complexvector adjoint_kernel(imsize);

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
                lensing_kernel[i * n + j] = std::complex<double>(real, imag);
                adjoint_kernel[i * n + j] = std::complex<double>(real, -imag);
            }
            else
            {
                lensing_kernel[i * n + j] = std::complex<double>(0, 0);
                adjoint_kernel[i * n + j] = std::complex<double>(0, 0);
            }
        }
    }

    auto forward = [=](complexvector &output, const complexvector &input)
    {
        for (int i = 0; i < (int)imsize; i++)
            output[i] = lensing_kernel[i] * input[i];
    };

    auto inverse = [=](complexvector &output, const complexvector &input)
    {
        output[0] = std::complex<double>(0., 0.);
        for (int i = 1; i < (int)imsize; i++)
            output[i] = input[i] / lensing_kernel[i];
    };

    auto adjoint = [=](complexvector &output, const complexvector &input)
    {
        for (int i = 0; i < (int)imsize; i++)
            output[i] = adjoint_kernel[i] * input[i];
    };

    return std::make_tuple(forward, inverse, adjoint);
}

void mmobservations::upsample(complexvector &hires, const complexvector &lowres)
{ // inputs and outputs in fourier space
    std::complex<double> _super(super, super);
    for (int i = 0; i < (int)(imsizey / 2); i++)
    {
        for (int j = 0; j < (int)(imsizex / 2); j++)
        {
            hires[i * superimsizey + j] = lowres[i * imsizey + j] * _super;

            hires[i * superimsizey + superimsizex - j - 1] = lowres[i * imsizey + imsizex - j - 1] * _super;

            hires[(superimsizey - i - 1) * superimsizey + j] = lowres[(imsizey - i - 1) * imsizey + j] * _super;

            hires[(superimsizey - i - 1) * superimsizey + superimsizex - j] = lowres[(imsizey - i - 1) * imsizey + imsizex - j - 1] * _super;
        }
    }
}

void mmobservations::downsample(complexvector &lowres, const complexvector &hires)
{ // inputs and outputs in fourier space
    std::complex<double> _super(super, super);
    for (int i = 0; i < (int)(imsizey / 2); i++)
    {
        for (int j = 0; j < (int)(imsizex / 2); j++)
        {
            lowres[i * imsizey + j] = hires[i * superimsizey + j] / _super;

            lowres[i * imsizey + imsizex - j - 1] = hires[i * superimsizey + superimsizex - j - 1] / _super;

            lowres[(imsizey - i - 1) * imsizey + j] = hires[(superimsizey - i - 1) * superimsizey + j] / _super;

            lowres[(imsizey - i - 1) * imsizey + imsizex - j] = hires[(superimsizey - i - 1) * superimsizey + superimsizex - j - 1] / _super;

        }
    }
}

void mmobservations::kaiser_squires(complexvector &output, const complexvector &input)
{
    complexvector temp(imsize);
    complexvector temp2(imsize);
    fft(temp, input);
    D(temp2, temp);
    ifft(output, temp2);
}

void mmobservations::kaiser_squires_inv(complexvector &output, const complexvector &input)
{
    complexvector temp(imsize);
    complexvector temp2(imsize);
    fft(temp, input);
    Dinv(temp2, temp);
    ifft(output, temp2);
}

void mmobservations::kaiser_squires_adj(complexvector &output, const complexvector &input)
{
    complexvector temp(imsize);
    complexvector temp2(imsize);
    fft(temp, input);
    Dadj(temp2, temp);
    ifft(output, temp2);
}
