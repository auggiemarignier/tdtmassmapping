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
                                   "tutorial_prior.txt",
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
    ASSERT_EQ(death.step(), 0);
    ASSERT_GT(death.propose, 0) << death.propose;
}

TEST_F(GlobalTest, Birth)
{
    BirthSliceMM birth(*global);
    ASSERT_EQ(birth.step(), 0);
    ASSERT_GT(birth.propose, 0) << birth.propose;
}

TEST_F(GlobalTest, Value)
{
    ValueSliceMM value(*global);
    ASSERT_EQ(value.step(), 0);
    ASSERT_GT(value.propose, 0) << value.propose;
}