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
    // Default constructor
    Observations(){};

#if 0 // Need to decide exactly how/when data gets read
    // Constructor that takes in vectors
    Observations(std::vector<double> _obs, std::vector<double> _sigma);

    // Constructor that takes a vector of obs and the stddev
    Observations(std::vector<double> _obs, double _sigma);

    // Constructor that takes a filename
    Observations(const char *filename);
#endif
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

class mmobservations : Observations
{
public:
    // Default constructor
    mmobservations(const uint _imsizex, const uint _imsizey);

#if 0 // Need to decide exactly how/when data gets read
    // Constructor that takes in vectors
    mmobservations(std::vector<double> _obs, std::vector<double> _sigma)
        : Observations(_obs, _sigma){};

    // Constructor that takes a vector of obs and the stddev
    mmobservations(std::vector<double> _obs, double _sigma)
        : Observations(_obs, _sigma){};

    // Constructor that takes a filename
    mmobservations(const char *filename)
        : Observations(filename){};
#endif

    double single_frequency_likelihood(
        complexvector model,
        std::complex<double> *residuals,
        std::complex<double> *residuals_normed,
        double &log_normalization);

    complexvector single_frequency_predictions(complexvector model);

    void kaiser_squires(fftw_complex *output, const fftw_complex *input);
    void kaiser_squires_inv(fftw_complex *output, const fftw_complex *input);
    void kaiser_squires_adj(fftw_complex *output, const fftw_complex *input);

    void set_observed_data(complexvector &_obs);
    void set_sigmas(std::vector<double> &_simgas);

private:
    std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>> init_fft_2d();

    std::function<void(fftw_complex *, const fftw_complex *)> fft;
    std::function<void(fftw_complex *, const fftw_complex *)> ifft;

    std::tuple<std::function<void(fftw_complex *, const fftw_complex *)>, std::function<void(fftw_complex *, const fftw_complex *)>,std::function<void(fftw_complex *, const fftw_complex *)>> build_lensing_kernels();

    std::function<void(fftw_complex *, const fftw_complex *)> D;
    std::function<void(fftw_complex *, const fftw_complex *)> Dadj;
    std::function<void(fftw_complex *, const fftw_complex *)> Dinv;

    const uint imsizex;
    const uint imsizey;
    const uint imsize;

    complexvector obs;
    std::vector<double> sigma;
};