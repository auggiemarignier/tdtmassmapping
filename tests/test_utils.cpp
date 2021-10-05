#include "gtest/gtest.h"

#include "utils.hpp"
#include "rng.hpp"

TEST(UtilsTest, MeanComplex)
{
    std::vector<std::complex<double>> vec(5);
    for (int i = 1; i <= vec.size(); i++)
        vec[i - 1] = std::complex<double>(i, vec.size() - i);
    auto mean = vector_mean(vec);
    ASSERT_FLOAT_EQ(mean.real(), 3.);
    ASSERT_FLOAT_EQ(mean.imag(), 2.);
};

TEST(UtilsTest, MeanDouble)
{
    std::vector<double> vec = {1., 2., 3., 4., 5.};
    auto mean = vector_mean(vec);
    ASSERT_FLOAT_EQ(mean, 3.);
};

TEST(UtilsTest, StdDevComplex)
{
    std::vector<std::complex<double>> vec(5);
    for (int i = 1; i <= vec.size(); i++)
        vec[i - 1] = std::complex<double>(i, vec.size() - i);
    auto stddev = vector_stddev(vec);
    ASSERT_FLOAT_EQ(stddev, 2.);
};

TEST(UtilsTest, StdDevDouble)
{
    std::vector<double> vec = {1., 2., 3., 4., 5.};
    auto stddev = vector_stddev(vec);
    ASSERT_FLOAT_EQ(stddev, std::sqrt(2.));
};

TEST(SatisticsTest, SNR)
{
    Rng random(1);
    std::vector<double> truth(50);
    std::vector<double> estimate(50);
    double l2_true = 0;
    for (int i = 0; i < truth.size(); i++)
    {
        truth[i] = random.normal(1.);
        estimate[i] = truth[i];
        l2_true += truth[i] * truth[i];
    }
    int rand_ind = random.uniform(truth.size());
    estimate[rand_ind] += 1;

    ASSERT_FLOAT_EQ(statistics::snr(truth, estimate), 10. * std::log10(l2_true));
};

TEST(StatisticsTest, Pearson)
{
    Rng random(1);
    std::vector<double> truth(50);
    std::vector<double> estimate(50);
    for (int i = 0; i < truth.size(); i++)
    {
        truth[i] = random.normal(1.);
        estimate[i] = truth[i];
    }

    ASSERT_FLOAT_EQ(statistics::pearson_correlation(truth, estimate), 1);
};

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}