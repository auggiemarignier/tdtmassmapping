#include "gtest/gtest.h"

#include "globalprop.hpp"
#include "mmobservations.hpp"
#include "proposals.hpp"

class BirthTest : public ::testing::Test
{
protected:
    mmobservations *observations;
    GlobalProposal *global;
    BirthProposal *birth;
    void SetUp() override
    {
        std::vector<double> obs = {0, 1, 2, 3};
        double sigma = 1.41;
        observations = new mmobservations(obs, sigma);

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
    ASSERT_EQ(global->mean_residual_n, 1);
}

class DeathTest : public ::testing::Test
{
protected:
    mmobservations *observations;
    GlobalProposal *global;
    DeathProposal *death;
    BirthProposal *birth;
    void SetUp() override
    {
        std::vector<double> obs = {0, 1, 2, 3};
        double sigma = 1.41;
        observations = new mmobservations(obs, sigma);

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
    ASSERT_EQ(global->mean_residual_n, 2);
}

class ValueTest : public ::testing::Test
{
protected:
    mmobservations *observations;
    GlobalProposal *global;
    ValueProposal *value;
    void SetUp() override
    {
        std::vector<double> obs = {0, 1, 2, 3};
        double sigma = 1.41;
        observations = new mmobservations(obs, sigma);

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
    ASSERT_EQ(global->mean_residual_n, 1);
}