#include "gtest/gtest.h"
#include <vector>
#include <array>
#include <cstring>
#include <cmath>
#include <algorithm>

#include "mmobservations.hpp"

class MMObsTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::vector<double> obs = {1, 2, 3, 4, 5};
        double sigma = 1.41;
        observations = new mmobservations(obs, sigma);
    }
    void TearDown() override
    {
        delete observations;
    }

    mmobservations *observations;
};

TEST_F(MMObsTest, MMPredsIdentity)
{
    std::vector<double> model = {1, 2, 3, 4, 5};
    ASSERT_EQ(model, observations->single_frequency_predictions(model));
    ASSERT_EQ(observations->n_obs, 5);
}

TEST_F(MMObsTest, MMLikelihood)
{
    std::vector<double> model = {1, 2, 3, 4, 5};
    independentgaussianhierarchicalmodel *hmodel = nullptr;
    hmodel = new independentgaussianhierarchicalmodel();
    hmodel->setparameter(0, 1.0);
    double log_normalization = 0.0;
    double residual[5];
    double residual_norm[5];

    ASSERT_EQ(0, observations->single_frequency_likelihood(
                     model, hmodel, residual, residual_norm, log_normalization));

    std::vector<double> model2 = {0, 2, 3, 4, 5};
    ASSERT_NE(0, observations->single_frequency_likelihood(
                     model2, hmodel, residual, residual_norm, log_normalization));
}

TEST_F(MMObsTest, FFTiFFT)
{
    const uint imsizex = 32;
    const uint imsizey = 32;
    auto fft_tuple = observations->init_fft_2d(imsizey, imsizex);
    auto fft = std::get<0>(fft_tuple);
    auto ifft = std::get<1>(fft_tuple);

    const double pi = std::acos(-1);
    const double freq = 2.;
    const double sps = 1. / 32.;

    std::array<std::complex<double>, imsizex * imsizey> f;
    for (uint j = 0; j < imsizex * imsizey; j++)
    {
        f[j] = std::complex<double>(std::sin(2 * pi * freq * sps * j), 0.);
    }

    fftw_complex *input = reinterpret_cast<fftw_complex *>(&f);
    fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsizex * imsizey);
    fftw_complex *recovered = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsizex * imsizey);

    fft(output, input);
    ifft(recovered, output);

    std::array<std::complex<double>, imsizex * imsizey> fhat;
    for (uint i = 1; i < imsizex * imsizey; i++)
    {
        // May fail at i=0 and imaginary component becuase 0 and 1e-30 are very different in ULP comparison 
        ASSERT_FLOAT_EQ(input[i][0], recovered[i][0]) << i;
        fhat[i] = std::complex<double>(output[i][0], output[i][1]);
    }

    uint index_max = 0;
    double current_max = std::norm(fhat[0]);
    for (uint i = 0; i < imsizex * imsizey; i++)
    {
        if (current_max < std::norm(fhat[i]))
        {
            current_max = std::norm(fhat[i]);
            index_max = i;
        }
    }
    ASSERT_EQ(index_max, (int)freq);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}