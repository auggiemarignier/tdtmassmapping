#include "gtest/gtest.h"
#include <vector>
#include <iostream>

#include "proposals.hpp"
#include "mmobservations.hpp"

class GlobalTest : public ::testing::Test
{
protected:
    GlobalSliceMM *global;
    void SetUp() override
    {
        global = new GlobalSliceMM("../../data/Bolshoi_7_clean_256.txt",
                                   NULL,
                                   8,
                                   8,
                                   1,
                                   100,
                                   4);
    }
};

TEST_F(GlobalTest, GlobalSetup)
{
    ASSERT_FLOAT_EQ(global->stddev, 0.034364982263599415);
}

TEST_F(GlobalTest, GlobalLikelihood)
{
    ASSERT_FALSE(true);
}