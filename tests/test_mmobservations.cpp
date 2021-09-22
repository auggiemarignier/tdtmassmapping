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

    // forward-backward recovery test
    std::array<std::complex<double>, imsize> fhat;
    for (uint i = 0; i < imsize; i++)
    {
        ASSERT_NEAR(input[i][0], recovered[i][0], 1e-12) << i;
        fhat[i] = std::complex<double>(output[i][0], output[i][1]);
    }

    // check spectrum peak is where it is expected
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

    // Parseval's Theorem to check unitary normalisation
    double space_power = 0.;
    double freq_power = 0.;
    for (uint i = 0; i < imsize; i ++)
    {
        space_power += std::norm(f[i]);
        freq_power += std::norm(fhat[i]);
    }
    ASSERT_NEAR(space_power, freq_power, 1e-12);
}

TEST_F(MMObsTest, LensingKernelInv)
{
    std::array<std::complex<double>, imsize> f;
    for (uint j = 0; j < imsize; j++)
    {
        f[j] = std::complex<double>(random->normal(1.), random->normal(1.));
    }

    fftw_complex *input = reinterpret_cast<fftw_complex *>(&f);
    fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *recovered = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->D(output, input);
    observations->Dinv(recovered, output);

    for (uint i = 1; i < imsize; i++)
    {
        // i = 0 can't be recovered due to mass sheet degeneracy
        ASSERT_NEAR(input[i][0], recovered[i][0], 1e-12) << i;
        ASSERT_NEAR(input[i][1], recovered[i][1], 1e-12) << i;
    }
}

TEST_F(MMObsTest, LensingKernelAdj)
{
    std::array<std::complex<double>, imsize> k1;
    std::array<std::complex<double>, imsize> g2;
    for (uint j = 0; j < imsize; j++)
    {
        k1[j] = std::complex<double>(random->normal(1.), random->normal(1.));
        g2[j] = std::complex<double>(random->normal(1.), random->normal(1.));
    }

    fftw_complex *kappa1 = reinterpret_cast<fftw_complex *>(&k1);
    fftw_complex *gamma2 = reinterpret_cast<fftw_complex *>(&g2);
    fftw_complex *kappa2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *gamma1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->D(gamma1, kappa1);
    observations->Dadj(kappa2, gamma2);

    std::array<std::complex<double>, imsize> k2;
    std::array<std::complex<double>, imsize> g1;
    for (uint j = 0; j < imsize; j++)
    {
        k2[j] = std::complex<double>(kappa2[j][0], kappa2[j][1]);
        g1[j] = std::complex<double>(gamma1[j][0], gamma1[j][1]);
    }

    std::complex<double> g1dotg2 = 0;
    std::complex<double> k1dotk2 = 0;
    for (uint j = 0; j < imsize; j++)
    {
        g1dotg2 += g1[j] * std::conj(g2[j]);
        k1dotk2 += k1[j] * std::conj(k2[j]);
    }

    ASSERT_FLOAT_EQ(g1dotg2.real(), k1dotk2.real());
    ASSERT_FLOAT_EQ(g1dotg2.imag(), k1dotk2.imag());
}

TEST_F(MMObsTest, KaiserSquiresInv)
{
    std::array<std::complex<double>, imsize> f;
    for (uint j = 0; j < imsize; j++)
    {
        f[j] = std::complex<double>(random->normal(1.), random->normal(1.));
    }

    fftw_complex *input = reinterpret_cast<fftw_complex *>(&f);
    fftw_complex *output = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *recovered = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->kaiser_squires(output, input);
    observations->kaiser_squires_inv(recovered, output);

    // Mass sheet degeneracy says we can only recover kappa up to a constant
    // Check input - recovered == constant
    double const_real = input[0][0] - recovered[0][0];
    double const_imag = input[0][1] - recovered[0][1];
    for (uint i = 1; i < imsize; i++)
    {
        ASSERT_NEAR(input[i][0] - recovered[i][0] - const_real, 0., 1e-12) << i;
        ASSERT_NEAR(input[i][1] - recovered[i][1] - const_imag, 0., 1e-12) << i;
    }
}

TEST_F(MMObsTest, KaiserSquiresAdj)
{
    std::array<std::complex<double>, imsize> k1;
    std::array<std::complex<double>, imsize> g2;
    for (uint j = 0; j < imsize; j++)
    {
        k1[j] = std::complex<double>(random->normal(1.), random->normal(1.));
        g2[j] = std::complex<double>(random->normal(1.), random->normal(1.));
    }

    fftw_complex *kappa1 = reinterpret_cast<fftw_complex *>(&k1);
    fftw_complex *gamma2 = reinterpret_cast<fftw_complex *>(&g2);
    fftw_complex *kappa2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *gamma1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->kaiser_squires(gamma1, kappa1);
    observations->kaiser_squires_adj(kappa2, gamma2);

    std::array<std::complex<double>, imsize> k2;
    std::array<std::complex<double>, imsize> g1;
    for (uint j = 0; j < imsize; j++)
    {
        k2[j] = std::complex<double>(kappa2[j][0], kappa2[j][1]);
        g1[j] = std::complex<double>(gamma1[j][0], gamma1[j][1]);
    }

    std::complex<double> g1dotg2 = 0;
    std::complex<double> k1dotk2 = 0;
    for (uint j = 0; j < imsize; j++)
    {
        g1dotg2 += g1[j] * std::conj(g2[j]);
        k1dotk2 += k1[j] * std::conj(k2[j]);
    }

    ASSERT_FLOAT_EQ(g1dotg2.real(), k1dotk2.real());
    ASSERT_FLOAT_EQ(g1dotg2.imag(), k1dotk2.imag());
}

TEST_F(MMObsTest, FFTAdj)
{
    std::array<std::complex<double>, imsize> k1;
    std::array<std::complex<double>, imsize> g2;
    for (uint j = 0; j < imsize; j++)
    {
        k1[j] = std::complex<double>(random->normal(1.), random->normal(1.));
        g2[j] = std::complex<double>(random->normal(1.), random->normal(1.));
    }

    fftw_complex *kappa1 = reinterpret_cast<fftw_complex *>(&k1);
    fftw_complex *gamma2 = reinterpret_cast<fftw_complex *>(&g2);
    fftw_complex *kappa2 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    fftw_complex *gamma1 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);

    observations->fft(gamma1, kappa1);
    observations->ifft(kappa2, gamma2);

    std::array<std::complex<double>, imsize> k2;
    std::array<std::complex<double>, imsize> g1;
    for (uint j = 0; j < imsize; j++)
    {
        k2[j] = std::complex<double>(kappa2[j][0], kappa2[j][1]);
        g1[j] = std::complex<double>(gamma1[j][0], gamma1[j][1]);
    }

    std::complex<double> g1dotg2 = 0;
    std::complex<double> k1dotk2 = 0;
    for (uint j = 0; j < imsize; j++)
    {
        g1dotg2 += g1[j] * std::conj(g2[j]);
        k1dotk2 += k1[j] * std::conj(k2[j]);
    }

    ASSERT_FLOAT_EQ(g1dotg2.real(), k1dotk2.real());
    ASSERT_FLOAT_EQ(g1dotg2.imag(), k1dotk2.imag());
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}