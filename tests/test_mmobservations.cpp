#include "gtest/gtest.h"
#include <vector>
#include <array>
#include <cstring>
#include <cmath>
#include <algorithm>

#include "mmobservations.hpp"
#include "rng.hpp"

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
    const uint imsize = imsizex * imsizey;
    auto fft_tuple = observations->init_fft_2d(imsizey, imsizex);
    auto fft = std::get<0>(fft_tuple);
    auto ifft = std::get<1>(fft_tuple);

    const double pi = std::acos(-1);
    const double freq = 2.;
    const double sps = 1. / 32.;

    std::array<std::complex<double>, imsize> f;
    for (uint j = 0; j < imsize; j++)
    {
        f[j] = std::complex<double>(std::sin(2 * pi * freq * sps * j), 0.);
    }

    fftw_complex *input = reinterpret_cast<fftw_complex *>(&f);
    fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *recovered = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    fft(output, input);
    ifft(recovered, output);

    std::array<std::complex<double>, imsize> fhat;
    for (uint i = 1; i < imsize; i++)
    {
        // May fail at i=0 and imaginary component becuase 0 and 1e-30 are very different in ULP comparison 
        ASSERT_FLOAT_EQ(input[i][0], recovered[i][0]) << i;
        fhat[i] = std::complex<double>(output[i][0], output[i][1]);
    }

    uint index_max = 0;
    double current_max = std::norm(fhat[0]);
    for (uint i = 1; i < imsize; i++)
    {
        if (current_max < std::norm(fhat[i]))
        {
            current_max = std::norm(fhat[i]);
            index_max = i;
        }
    }
    ASSERT_EQ(index_max, (int)freq);
}

TEST_F(MMObsTest, LensingKernel)
{
    const uint imsizex = 32;
    const uint imsizey = 32;
    const uint imsize = imsizex * imsizey;
    auto operator_tuple = observations->build_lensing_kernels(imsizex, imsizey);
    auto D = std::get<0>(operator_tuple);
    auto D_adj = std::get<1>(operator_tuple);

    Rng random(1);

    std::array<std::complex<double>, imsize> f;
    for (uint j = 0; j < imsize; j++)
    {
        f[j] = std::complex<double>(random.uniform(), random.uniform());
    }

    fftw_complex *input = reinterpret_cast<fftw_complex *>(&f);
    fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *recovered = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    D(output, input);
    D_adj(recovered, output);

    for (uint i = 0; i < imsize; i++)
    {
        ASSERT_FLOAT_EQ(input[i][0], recovered[i][0]) << i;
        ASSERT_FLOAT_EQ(input[i][1], recovered[i][1]) << i;
        
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}