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
        std::vector<double> obs = {1, 2, 3, 4, 5};
        std::vector<double> sigma = {1.41, 1.41, 1.41, 1.41, 1.41};
        global = new GlobalSliceMM(obs,
                                   sigma,
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
    ASSERT_FLOAT_EQ(global->observations->sigma[0], 1.41);
}

TEST_F(GlobalTest, GlobalLikelihood)
{
    // wavetree model at initialisation is 0 everywhere
    double likelihood = global->likelihood(global->current_log_normalization);
    ASSERT_FLOAT_EQ(likelihood, 13.8323);
}