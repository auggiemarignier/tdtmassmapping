#include "gtest/gtest.h"

#include "globalprop.hpp"
#include "mmobservations.hpp"
#include "proposals.hpp"

class ProposalTest : public ::testing::Test
{
protected:
    mmobservations *observations;
    GlobalProposal *global;
    Proposal *proposal;
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

        proposal = new Proposal(*global, "Base");
    }
};