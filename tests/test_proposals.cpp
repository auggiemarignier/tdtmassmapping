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
    mmobservations *observations;
    GlobalProposal *global;
    DeathProposal *death;
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

        death = new DeathProposal(*global);
    }
};

TEST_F(DeathTest, DeathStep)
{
    ASSERT_EQ(death->name, 2);
    int step_taken = death->step();
    ASSERT_EQ(death->propose, 1);
    ASSERT_TRUE(step_taken == 0 || step_taken == 1);
}
