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

TEST_F(GlobalTest, GlobalSetup)
{
    ASSERT_FLOAT_EQ(global->stddev, 0.034364982263599415);
}
