#include "gtest/gtest.h"
#include <vector>

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
    mmobservations observations;
};

TEST_F(MMObsTest, MMPredsIdentity)
{
    std::vector<double> model = {1, 2, 3, 4, 5};
    ASSERT_EQ(model, observations.single_frequency_predictions(model));
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}