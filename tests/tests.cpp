#include "gtest/gtest.h"

#include "proposals.hpp"
#include "mmobservations.hpp"

class GlobalTest : public ::testing::Test
{
protected:
    GlobalSliceMM *global;
    void SetUp() override
    {
        global = new GlobalSliceMM("../../data/Bolshoi_7_clean_256.txt");
    }

    void TearDown() override
    {
        delete global;
    }
};

class MMObsTest : public ::testing::Test
{
protected:
    mmobservations *observations;
};

TEST_F(MMObsTest, MMPredsIdentity)
{
    ASSERT_FALSE(observations->single_frequency_predictions(nullptr));

    const double model1 = 1.0;
    ASSERT_TRUE(observations->single_frequency_predictions(&model1));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}