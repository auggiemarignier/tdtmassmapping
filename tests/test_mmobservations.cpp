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

TEST_F(MMObsTest, LensingKernelInv)
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
    observations->Dinv(recovered, output);

    for (uint i = 1; i < imsize; i++)
    {
        // i = 0 can't be recovered due to mass sheet degeneracy
        ASSERT_FLOAT_EQ(input[i][0], recovered[i][0]) << i;
        ASSERT_FLOAT_EQ(input[i][1], recovered[i][1]) << i;
    }
}

TEST_F(MMObsTest, LensingKernelAdj)
{
    std::array<std::complex<double>, imsize> k1;
    std::array<std::complex<double>, imsize> k2;
    for (uint j = 0; j < imsize; j++)
    {
        k1[j] = std::complex<double>(random->uniform(), random->uniform());
        k2[j] = std::complex<double>(random->uniform(), random->uniform());
    }

    fftw_complex *kappa1 = reinterpret_cast<fftw_complex *>(&k1);
    fftw_complex *kappa2 = reinterpret_cast<fftw_complex *>(&k2);
    fftw_complex *gamma1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *gamma2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *kappa2prime = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->D(gamma1, kappa1);
    observations->D(gamma2, kappa2);
    observations->Dadj(kappa2prime, gamma2);

    std::array<std::complex<double>, imsize> g1;
    std::array<std::complex<double>, imsize> g2;
    std::array<std::complex<double>, imsize> k2p;
    for (uint j = 0; j < imsize; j++)
    {
        g1[j] = std::complex<double>(gamma1[j][0], gamma1[j][1]);
        g2[j] = std::complex<double>(gamma2[j][0], gamma2[j][1]);
        k2p[j] = std::complex<double>(kappa2prime[j][0], kappa2prime[j][1]);
    }

    double g1dotg2 = 0;
    double k1dotk2p = 0;
    for (uint j = 0; j < imsize; j++)
    {
        g1dotg2 += std::abs(g1[j] * g2[j]);
        k1dotk2p += std::abs(k1[j] * k2p[j]);
    }

    ASSERT_FLOAT_EQ(g1dotg2, k1dotk2p);
}

TEST_F(MMObsTest, KaiserSquiresInv)
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
    observations->kaiser_squires_inv(recovered, output);

    for (uint i = 1; i < imsize; i++)
    {
        // i = 0 can't be recovered due to mass sheet degeneracy
        ASSERT_FLOAT_EQ(input[i][0], recovered[i][0]) << i;
        ASSERT_FLOAT_EQ(input[i][1], recovered[i][1]) << i;
    }
}

TEST_F(MMObsTest, KaiserSquiresAdj)
{
    std::array<std::complex<double>, imsize> k1;
    std::array<std::complex<double>, imsize> k2;
    for (uint j = 0; j < imsize; j++)
    {
        k1[j] = std::complex<double>(random->uniform(), random->uniform());
        k2[j] = std::complex<double>(random->uniform(), random->uniform());
    }

    fftw_complex *kappa1 = reinterpret_cast<fftw_complex *>(&k1);
    fftw_complex *kappa2 = reinterpret_cast<fftw_complex *>(&k2);
    fftw_complex *gamma1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *gamma2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *kappa2prime = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->kaiser_squires(gamma1, kappa1);
    observations->kaiser_squires(gamma2, kappa2);
    observations->kaiser_squires_adj(kappa2prime, gamma2);

    std::array<std::complex<double>, imsize> g1;
    std::array<std::complex<double>, imsize> g2;
    std::array<std::complex<double>, imsize> k2p;
    for (uint j = 0; j < imsize; j++)
    {
        g1[j] = std::complex<double>(gamma1[j][0], gamma1[j][1]);
        g2[j] = std::complex<double>(gamma2[j][0], gamma2[j][1]);
        k2p[j] = std::complex<double>(kappa2prime[j][0], kappa2prime[j][1]);
    }

    double g1dotg2 = 0;
    double k1dotk2p = 0;
    for (uint j = 0; j < imsize; j++)
    {
        g1dotg2 += std::abs(g1[j] * g2[j]);
        k1dotk2p += std::abs(k1[j] * k2p[j]);
    }

    ASSERT_FLOAT_EQ(g1dotg2, k1dotk2p);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}