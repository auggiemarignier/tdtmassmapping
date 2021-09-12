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
        const uint imsizex = 32;
        const uint imsizey = 32;
        observations = new mmobservations(imsizex, imsizey);
        random = new Rng(1);
    }
    void TearDown() override
    {
        delete observations;
    }

    mmobservations *observations;
    Rng *random;
    static constexpr uint imsize = 32 * 32;
};

TEST_F(MMObsTest, FFTiFFT)
{
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

    observations->fft(output, input);
    observations->ifft(recovered, output);

    std::array<std::complex<double>, imsize> fhat;
    for (uint i = 1; i < imsize; i++)
    {
        // May fail at i=0 and imaginary component becuase 0 and 1e-30 are very different in ULP comparison
        ASSERT_FLOAT_EQ(input[i][0], recovered[i][0]) << i;
        fhat[i] = std::complex<double>(output[i][0], output[i][1]);
    }

    int index_max = 0;
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
    std::array<std::complex<double>, imsize> f;
    for (uint j = 0; j < imsize; j++)
    {
        f[j] = std::complex<double>(random->uniform(), random->uniform());
    }

    fftw_complex *input = reinterpret_cast<fftw_complex *>(&f);
    fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *recovered = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->D(output, input);
    observations->Dadj(recovered, output);

    for (uint i = 0; i < imsize; i++)
    {
        ASSERT_FLOAT_EQ(input[i][0], recovered[i][0]) << i;
        ASSERT_FLOAT_EQ(input[i][1], recovered[i][1]) << i;
    }
}

TEST_F(MMObsTest, KaiserSquires)
{
    std::array<std::complex<double>, imsize> f;
    for (uint j = 0; j < imsize; j++)
    {
        f[j] = std::complex<double>(random->uniform(), random->uniform());
    }

    fftw_complex *input = reinterpret_cast<fftw_complex *>(&f);
    fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *recovered = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->kaiser_squires(output, input);
    observations->kaiser_squires_adj(recovered, output);

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