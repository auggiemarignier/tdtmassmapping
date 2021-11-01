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
        observations = new mmobservations(imsizex, imsizey, 1);
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

TEST(MMObsDummyTest, UpsampleDownsample)
{
    const uint imsizex = 32;
    const uint imsizey = 32;
    const uint imsize = imsizex * imsizey;
    const uint super = 2;
    const uint superimsizex = super * imsizex;
    const uint superimsizey = super * imsizey;
    const uint superimsize = super * imsize;
    const double pi = 4. * acos(1. / sqrt(2));

    mmobservations observations(imsizex, imsizey, super);
    complexvector kappa(superimsize);
    for (int i = 0; i < superimsizey; i++)
    {
        for (int j = 0; j < superimsizex; j++)
        {
            kappa[i * superimsizey + j] = std::complex<double>(
                sin(i / pi) + sin(j / pi),
                sin(2 * i / pi) + sin(2 * j / pi));
            kappa[i * superimsizey + j] += std::complex<double>(
                sin(i / (2 * pi)) + sin(j / (2 * pi)),
                -sin(i / pi) - sin(j / pi));
        }
    }

    FILE *fp = fopen("/Users/auggiemarignier/Documents/PhD/TDT/massmapping/outputs/checker0.txt", "w");
    for (int i = 0; i < superimsizey; i++)
    {
        for (int j = 0; j < superimsizex; j++)
        {
            fprintf(fp, "%10.6f ", kappa[i * superimsizey + j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    complexvector kappahat(superimsize);
    complexvector kappadownhat(imsize);
    complexvector kappadown(imsize);
    complexvector kapparechat(superimsize);
    complexvector kapparec(superimsize);

    observations.s_fft(kappahat, kappa);
    observations.downsample(kappadownhat, kappahat);
    observations.ifft(kappadown, kappadownhat);
    observations.upsample(kapparechat, kappadownhat);
    observations.s_ifft(kapparec, kapparechat);

    FILE *fp0 = fopen("/Users/auggiemarignier/Documents/PhD/TDT/massmapping/outputs/checker.txt", "w");
    FILE *fp1 = fopen("/Users/auggiemarignier/Documents/PhD/TDT/massmapping/outputs/checkerdown.txt", "w");
    FILE *fp2 = fopen("/Users/auggiemarignier/Documents/PhD/TDT/massmapping/outputs/checkerrec.txt", "w");

    for (int i = 0; i < imsizey; i++)
    {
        for (int j = 0; j < imsizex; j++)
        {
            fprintf(fp1, "%10.6f ", kappadown[i * imsizey + j].real());
        }
        fprintf(fp1, "\n");
    }
    for (int i = 0; i < (int)(superimsizey); i++)
    {
        for (int j = 0; j < (int)(superimsizex); j++)
        {
            fprintf(fp0, "%10.6f ", kappa[i * superimsizey + j].real());
            fprintf(fp2, "%10.6f ", kapparec[i * superimsizey + j].real());
        }
        fprintf(fp0, "\n");
        fprintf(fp2, "\n");
    }

    fclose(fp0);
    fclose(fp1);
    fclose(fp2);

    for (uint i = 0; i < imsize; i++)
    {
        EXPECT_NEAR(kappa[i].real() - kapparec[i].real(), 0., 1e-12) << i;
        EXPECT_NEAR(kappa[i].imag() - kapparec[i].imag(), 0., 1e-12) << i;
    }
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}