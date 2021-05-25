#include "gtest/gtest.h"
#include <vector>
#include <iostream>

#include "proposals.hpp"
#include "mmobservations.hpp"

#include "wavetomo2dutil.hpp"

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

TEST_F(GlobalTest, GlobalAccept)
{
    ASSERT_EQ(global->mean_residual_n, 0);
    for (int i = 0; i < global->residual_size; i++)
    {
        ASSERT_EQ(global->mean_residual[i], 0);
        ASSERT_EQ(global->last_valid_residual[i], 0);
    }

    double likelihood = global->likelihood(global->current_log_normalization);
    global->accept();

    ASSERT_EQ(global->mean_residual_n, 1);
    for (int i = 0; i < global->residual_size; i++)
    {
        ASSERT_EQ(global->mean_residual[i], global->observations->obs[i]);
        ASSERT_EQ(global->last_valid_residual[i], global->observations->obs[i]);
    }
}

TEST_F(GlobalTest, Death)
{
    DeathSliceMM death(*global);
    for (int i = 0; i < global->treemaxdepth; i++)
    {
        ASSERT_EQ(death.propose_depth[i], 0) << i << " " << death.propose_depth[i];
    }
    ASSERT_EQ(death.propose, 0) << death.propose;
    ASSERT_EQ(death.step(), 0) << death.step();
    ASSERT_GT(death.propose, 0) << death.propose;
    for (int i = 0; i < global->treemaxdepth; i++)
    {
        ASSERT_EQ(death.propose_depth[i], 0) << i << " " << death.propose_depth[i];
    }
    ASSERT_FALSE(true) << death.global.treemaxdepth;
}