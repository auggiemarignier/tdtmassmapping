#include "gtest/gtest.h"

#include "globalprop.hpp"
#include "mmobservations.hpp"
#include "proposals.hpp"

class BirthTest : public ::testing::Test
{
protected:
    Identity *observations;
    GlobalProposal *global;
    BirthProposal *birth;
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
        observations->set_data_errors(sigma);

        global = new GlobalProposal(observations,
                                    NULL,
                                    "tutorial_prior.txt",
                                    1,
                                    1,
                                    1,
                                    4,
                                    4);

        birth = new BirthProposal(*global);
    }
};

TEST_F(BirthTest, BirthStep)
{
    ASSERT_EQ(birth->name, 1);
    int step_taken = birth->step();
    ASSERT_EQ(birth->propose, 1);
    ASSERT_TRUE(step_taken == 0 || step_taken == 1);
}

class DeathTest : public ::testing::Test
{
protected:
    Identity *observations;
    GlobalProposal *global;
    DeathProposal *death;
    BirthProposal *birth;
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
        observations->set_data_errors(sigma);

        global = new GlobalProposal(observations,
                                    NULL,
                                    "tutorial_prior.txt",
                                    1,
                                    1,
                                    1,
                                    4,
                                    4);

        death = new DeathProposal(*global);

        // need to birth before you can kill a node
        birth = new BirthProposal(*global);
        birth->step();
    }
};

TEST_F(DeathTest, DeathStep)
{
    ASSERT_EQ(death->name, 2);
    int step_taken = death->step();
    ASSERT_EQ(death->propose, 1);
    ASSERT_TRUE(step_taken == 0 || step_taken == 1);
}

class ValueTest : public ::testing::Test
{
protected:
    Identity *observations;
    GlobalProposal *global;
    ValueProposal *value;
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
        observations->set_data_errors(sigma);

        global = new GlobalProposal(observations,
                                    NULL,
                                    "tutorial_prior.txt",
                                    1,
                                    1,
                                    1,
                                    4,
                                    4);

        value = new ValueProposal(*global);
    }
};

TEST_F(ValueTest, ValueStep)
{
    ASSERT_EQ(value->name, 3);
    int step_taken = value->step();
    ASSERT_EQ(value->propose, 1);
    ASSERT_TRUE(step_taken == 0 || step_taken == 1);
}