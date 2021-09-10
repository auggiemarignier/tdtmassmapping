#pragma once

#include <stdio.h>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <tuple>
#include <cassert>
#include <functional>
#include <memory>

#include "hierarchicalmodel.hpp"

typedef std::vector<std::complex<double>> complexvector;
typedef unsigned int uint;

class Observations
{
public:
    // Constructor that takes in vectors
    Observations(std::vector<double> _obs, std::vector<double> _sigma);

    // Constructor that takes a vector of obs and the stddev
    Observations(std::vector<double> _obs, double _sigma);

    // Constructor that takes a filename
    Observations(const char *filename);

    virtual ~Observations(){};

    virtual double single_frequency_likelihood(std::vector<double> model,
                                               const hierarchicalmodel *hmodel,
                                               double *residuals,
                                               double *residuals_normed,
                                               double &log_normalization) = 0;

    virtual std::vector<double> single_frequency_predictions(std::vector<double> model) = 0;

    virtual bool save_residuals(const char *filename,
                                const double *residuals,
                                const double *residuals_normed);

    std::vector<double> obs;
    std::vector<double> sigma;
    size_t n_obs;
};

class mmobservations : public Observations
{
public:
    // Constructor that takes in vectors
    mmobservations(std::vector<double> _obs, std::vector<double> _sigma)
        : Observations(_obs, _sigma){};

    // Constructor that takes a vector of obs and the stddev
    mmobservations(std::vector<double> _obs, double _sigma)
        : Observations(_obs, _sigma){};

    // Constructor that takes a filename
    mmobservations(const char *filename)
        : Observations(filename){};

    double single_frequency_likelihood(
        std::vector<double> model,
        const hierarchicalmodel *hmodel,
        double *residuals,
        double *residuals_normed,
        double &log_normalization) override;

    std::vector<double> single_frequency_predictions(std::vector<double> model) override;

    std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>> init_fft_2d(const uint &imsizey, const uint &imsizex);

    std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>> build_lensing_kernels(const uint &imsizey, const uint &imsizex);

private:
    std::shared_ptr<fftw_plan_s> plan_forward;
    std::shared_ptr<fftw_plan_s> plan_inverse;

    complexvector lensing_kernel, adjoint_kernel;

};