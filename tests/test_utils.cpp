#include "gtest/gtest.h"

#include "utils.hpp"

TEST(UtilsTest, TestMean)
{
    std::vector<std::complex<double>> vec(5);
    for (int i = 1; i <= vec.size(); i++)
        vec[i-1] = std::complex<double>(i, vec.size() - i);
    auto mean = vector_mean(vec);
    ASSERT_FLOAT_EQ(mean.real(), 3.);
    ASSERT_FLOAT_EQ(mean.imag(), 2.);
};

TEST(UtilsTest, TestStdDev)
{
    std::vector<std::complex<double>> vec(5);
    for (int i = 1; i <= vec.size(); i++)
        vec[i-1] = std::complex<double>(i, vec.size() - i);
    auto stddev = vector_stddev(vec);
    ASSERT_FLOAT_EQ(stddev, 2.);
};

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}