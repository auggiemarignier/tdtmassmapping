#include "gtest/gtest.h"

#include "utils.hpp"

TEST(UtilsTest, TestMeanComplex)
{
    std::vector<std::complex<double>> vec(5);
    for (int i = 1; i <= vec.size(); i++)
        vec[i-1] = std::complex<double>(i, vec.size() - i);
    auto mean = vector_mean(vec);
    ASSERT_FLOAT_EQ(mean.real(), 3.);
    ASSERT_FLOAT_EQ(mean.imag(), 2.);
};

TEST(UtilsTest, TestMeanDouble)
{
    std::vector<double> vec = { 1., 2., 3., 4., 5.};
    auto mean = vector_mean(vec);
    ASSERT_FLOAT_EQ(mean, 3.);
};

TEST(UtilsTest, TestStdDevComplex)
{
    std::vector<std::complex<double>> vec(5);
    for (int i = 1; i <= vec.size(); i++)
        vec[i-1] = std::complex<double>(i, vec.size() - i);
    auto stddev = vector_stddev(vec);
    ASSERT_FLOAT_EQ(stddev, 2.);
};

TEST(UtilsTest, TestStdDevDouble)
{
    std::vector<double> vec = { 1., 2., 3., 4., 5.};
    auto stddev = vector_stddev(vec);
    ASSERT_FLOAT_EQ(stddev, std::sqrt(2.));
};

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}