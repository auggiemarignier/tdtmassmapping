#include "gtest/gtest.h"

#include "globalprop.hpp"
#include "mmobservations.hpp"

class GlobalPropTest : public ::testing::Test
{
protected:
    mmobservations *observations;
    GlobalProposal *global;
    void SetUp() override
    {
        std::vector<double> obs = {1, 2, 3, 4, 5};
        double sigma = 1.41;
        observations = new mmobservations(obs, sigma);

        global = new GlobalProposal(NULL,
                                    observations,
                                    NULL,
                                    "tutorial_prior.txt",
                                    8,
                                    8,
                                    1,
                                    100,
                                    4);
    }
};

TEST_F(GlobalPropTest, GlobalSetup)
{
    ASSERT_FLOAT_EQ(global->observations->sigma[0], 1.41);
}

TEST_F(GlobalPropTest, GlobalLikelihood)
{
    // wavetree model at initialisation is 0 everywhere
    double likelihood = global->likelihood(global->current_log_normalization);
    ASSERT_FLOAT_EQ(likelihood, 13.8323);
}

TEST_F(GlobalPropTest, GlobalAccept)
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