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
        complexvector obs = {std::complex<double>(1., 0.),
                             std::complex<double>(2., 0.),
                             std::complex<double>(3., 0.),
                             std::complex<double>(4., 0.),
                             std::complex<double>(5., 0.)};
        std::vector<double> sigma = {1.41};
        observations = new Identity();
        observations->set_observed_data(obs);
        observations->set_data_errors(sigma, false);

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

TEST_F(GlobalPropTest, GlobalPrior)
{
    // wavetree model at initialisation is 0 everywhere with k=1
    // only 1 possible tree arrangement arrangement when k=1
    // initial tree values given 0 prior probability
    double p = global->prior();
    double expected = log(1. / global->kmax) + log(1.) + 0;
    ASSERT_EQ(p, expected);
}