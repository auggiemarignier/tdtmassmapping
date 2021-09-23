#include "gtest/gtest.h"

#include "globalprop.hpp"
#include "mmobservations.hpp"

class GlobalPropTest : public ::testing::Test
{
protected:
    Identity *observations;
    GlobalProposal *global;
    void SetUp() override
    {
        std::vector<double> obs = {1, 2, 3, 4, 5};
        double sigma = 1.41;
        observations = new Identity();

        global = new GlobalProposal(observations,
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
