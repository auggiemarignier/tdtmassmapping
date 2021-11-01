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
        const int super = 1;
        observations = new mmobservations(imsizex, imsizey, super);
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

TEST_F(MMObsTest, MMLikelihood)
{
    complexvector f(imsize);
    ;
    for (uint j = 0; j < imsize; j++)
        f[j] = std::complex<double>(random->normal(1.), 0);

    complexvector predictions = observations->single_frequency_predictions(f);
    for (int i = 1; i < imsize; i++)
        ASSERT_NE(std::complex<double>(0, 0), predictions[i]);
    observations->set_observed_data(predictions);

    std::vector<double> sig = {1.};
    observations->set_sigmas(sig);

    double log_normalization = 0.0;

    ASSERT_EQ(0, observations->single_frequency_likelihood(f, log_normalization));
}

TEST_F(MMObsTest, KaiserSquiresInv)
{
    complexvector input(imsize);
    complexvector output(imsize);
    complexvector recovered(imsize);
    for (uint j = 0; j < imsize; j++)
        input[j] = std::complex<double>(random->normal(1.), random->normal(1.));

    observations->kaiser_squires(output, input);
    observations->kaiser_squires_inv(recovered, output);
    for (int i = 1; i < imsize; i++)
    {
        ASSERT_NE(std::complex<double>(0, 0), output[i]);
        ASSERT_NE(std::complex<double>(0, 0), recovered[i]);
    }
    // Mass sheet degeneracy says we can only recover kappa up to a constant
    // Check input - recovered == constant
    std::complex<double> const_diff = input[0] - recovered[0];
    for (uint i = 1; i < imsize; i++)
    {
        std::complex<double> diff = input[i] - recovered[i];
        ASSERT_NEAR(diff.real() - const_diff.real(), 0., 1e-12) << i;
        ASSERT_NEAR(diff.imag() - const_diff.imag(), 0., 1e-12) << i;
    }
}

TEST_F(MMObsTest, KaiserSquiresAdj)
{
    complexvector k1(imsize);
    complexvector k2(imsize);
    complexvector g1(imsize);
    complexvector g2(imsize);
    for (uint j = 0; j < imsize; j++)
    {
        k1[j] = std::complex<double>(random->normal(1.), random->normal(1.));
        g2[j] = std::complex<double>(random->normal(1.), random->normal(1.));
    }

    observations->kaiser_squires(g1, k1);
    observations->kaiser_squires_adj(k2, g2);
    for (int i = 1; i < imsize; i++)
    {
        ASSERT_NE(std::complex<double>(0, 0), g1[i]);
        ASSERT_NE(std::complex<double>(0, 0), k2[i]);
    }

    std::complex<double> g1dotg2 = 0;
    std::complex<double> k1dotk2 = 0;
    for (uint j = 0; j < imsize; j++)
    {
        g1dotg2 += g1[j] * std::conj(g2[j]);
        k1dotk2 += k1[j] * std::conj(k2[j]);
    }
    ASSERT_NE(std::complex<double>(0, 0), g1dotg2);
    ASSERT_NE(std::complex<double>(0, 0), k1dotk2);
    ASSERT_FLOAT_EQ(g1dotg2.real(), k1dotk2.real());
    ASSERT_FLOAT_EQ(g1dotg2.imag(), k1dotk2.imag());
}

class MMObsSuperTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        const uint imsizex = 32;
        const uint imsizey = 32;
        const int super = 2;
        observations = new mmobservations(imsizex, imsizey, super);
        random = new Rng(1);
    }
    void TearDown() override
    {
        delete observations;
    }

    mmobservations *observations;
    Rng *random;
    static constexpr uint imsize = 32 * 32;
    static constexpr uint superimsize = 64 * 64;
};

TEST_F(MMObsSuperTest, KaiserSquiresInv)
{
    complexvector input(superimsize);
    complexvector output(imsize);
    complexvector recovered(superimsize);
    for (uint j = 0; j < superimsize; j++)
        input[j] = std::complex<double>(random->normal(1.), random->normal(1.));

    observations->kaiser_squires(output, input);
    observations->kaiser_squires_inv(recovered, output);
    for (int i = 1; i < imsize; i++)
        ASSERT_NE(std::complex<double>(0, 0), output[i]);
    for (int i = 1; i < superimsize; i++)
        ASSERT_NE(std::complex<double>(0, 0), recovered[i]);

    // Mass sheet degeneracy says we can only recover kappa up to a constant
    // Check input - recovered == constant
    std::complex<double> const_diff = input[0] - recovered[0];
    for (uint i = 1; i < superimsize; i++)
    {
        std::complex<double> diff = input[i] - recovered[i];
        ASSERT_NEAR(diff.real() - const_diff.real(), 0., 1e-12) << i;
        ASSERT_NEAR(diff.imag() - const_diff.imag(), 0., 1e-12) << i;
    }
}

TEST_F(MMObsSuperTest, KaiserSquiresAdj)
{
    complexvector k1(superimsize);
    complexvector k2(superimsize);
    complexvector g1(imsize);
    complexvector g2(imsize);
    for (uint j = 0; j < superimsize; j++)
        k1[j] = std::complex<double>(random->normal(1.), random->normal(1.));
    for (uint j = 0; j < imsize; j++)
        g2[j] = std::complex<double>(random->normal(1.), random->normal(1.));

    observations->kaiser_squires(g1, k1);
    observations->kaiser_squires_adj(k2, g2);
    for (int i = 1; i < imsize; i++)
        ASSERT_NE(std::complex<double>(0, 0), g1[i]);
    for (int i = 1; i < superimsize; i++)
        ASSERT_NE(std::complex<double>(0, 0), k2[i]);

    std::complex<double> g1dotg2 = 0;
    std::complex<double> k1dotk2 = 0;
    for (uint j = 0; j < imsize; j++)
        g1dotg2 += g1[j] * std::conj(g2[j]);
    for (uint j = 0; j < superimsize; j++)
        k1dotk2 += k1[j] * std::conj(k2[j]);

    ASSERT_NE(std::complex<double>(0, 0), g1dotg2);
    ASSERT_NE(std::complex<double>(0, 0), k1dotk2);
    ASSERT_FLOAT_EQ(g1dotg2.real(), k1dotk2.real());
    ASSERT_FLOAT_EQ(g1dotg2.imag(), k1dotk2.imag());
}
int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}