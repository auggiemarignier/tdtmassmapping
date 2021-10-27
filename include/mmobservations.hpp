#pragma once

#include <stdio.h>
#include <vector>
#include <complex>
#include <fftw3.h>
#include <tuple>
#include <cassert>
#include <functional>
#include <memory>

typedef std::vector<std::complex<double>> complexvector;
typedef unsigned int uint;

class Observations
{
public:
    // Default constructor
    Observations(){};

#if 0 // Need to decide exactly how/when data gets read
    // Constructor that takes in vectors
    Observations(complexvector _obs, complexvector _sigma);

    // Constructor that takes a vector of obs and the stddev
    Observations(complexvector _obs, double _sigma);

    // Constructor that takes a filename
    Observations(const char *filename);
#endif
    virtual ~Observations(){};

    virtual double single_frequency_likelihood(complexvector &model,
                                               double &log_normalization);

    virtual complexvector single_frequency_predictions(complexvector &model) = 0;

    virtual void set_observed_data(complexvector &_obs);
    virtual void set_sigmas(std::vector<double> &_simgas);

    complexvector obs;
    std::vector<double> sigma;
    size_t n_obs;
};

class Identity : public Observations
{
public:
    complexvector single_frequency_predictions(complexvector &model) override;
};

class mmobservations : public Observations
{
public:
    // Default constructor
    mmobservations(const uint _imsizex, const uint _imsizey, const uint super);

#if 0 // Need to decide exactly how/when data gets read
    // Constructor that takes in vectors
    mmobservations(complexvector _obs, complexvector _sigma)
        : Observations(_obs, _sigma){};

    // Constructor that takes a vector of obs and the stddev
    mmobservations(complexvector _obs, double _sigma)
        : Observations(_obs, _sigma){};

    // Constructor that takes a filename
    mmobservations(const char *filename)
        : Observations(filename){};
#endif

    complexvector single_frequency_predictions(complexvector &model) override;

    void kaiser_squires(complexvector &output, const complexvector &input);
    void kaiser_squires_inv(complexvector &output, const complexvector &input);
    void kaiser_squires_adj(complexvector &output, const complexvector &input);

private:
    std::tuple<std::function<void(complexvector &, const complexvector &)>, std::function<void(complexvector &, const complexvector &)>> init_fft_2d(const uint _imsizex, const uint _imsizey);

    std::function<void(complexvector &, const complexvector &)> fft;
    std::function<void(complexvector &, const complexvector &)> ifft;
    std::function<void(complexvector &, const complexvector &)> s_fft;
    std::function<void(complexvector &, const complexvector &)> s_ifft;

    std::tuple<std::function<void(complexvector &, const complexvector &)>, std::function<void(complexvector &, const complexvector &)>, std::function<void(complexvector &, const complexvector &)>> build_lensing_kernels();

    std::function<void(complexvector &, const complexvector &)> D;
    std::function<void(complexvector &, const complexvector &)> Dadj;
    std::function<void(complexvector &, const complexvector &)> Dinv;

    const uint imsizex;
    const uint imsizey;
    const uint imsize;
    const uint super;
    uint superimsizex;
    uint superimsizey;
    uint superimsize;
};