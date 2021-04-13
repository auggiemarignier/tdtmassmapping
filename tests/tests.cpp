#include "gtest/gtest.h"

#include "mmobservations.hpp"

TEST(MMObsTest, MMPredsIdentity)
{
    const double model = 0.0;
    ASSERT_FALSE(model);

    const double model1 = 1.0;
    ASSERT_TRUE(model1);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}