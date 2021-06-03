#pragma once

#include "globalprop.hpp"

#include <mpi.h>

class Proposals
{
public:
    Proposal(GlobalProposal &global);
    virtual ~Proposal()
    {
        delete[] propose_depth;
        delete[] accept_depth;
    };

    virtual int step() = 0;

    std::string write_short_stats();

    std::string write_long_stats();

    void initialize_mpi(MPI_Comm communicator);

    GlobalProposal &global;
    int propose;
    int accept;

    int *propose_depth;
    int *accept_depth;

private:
    virtual int choose_proposal_location_and_value(int k,
                                                   double &ratio,
                                                   int &prop_depth,
                                                   int &prop_idx,
                                                   double &choose_prob,
                                                   double &prop_value,
                                                   double &prop_prob,
                                                   int &prop_valid,
                                                   double &prop_parent_coeff,
                                                   int &ii,
                                                   int &ij);

    virtual int communicate_proposal_location_and_value(int &prop_valid,
                                                        int &prop_idx,
                                                        int &prop_depth,
                                                        double &prop_value);

    virtual int propose_proposal(int &prop_valid,
                                 int &prop_idx,
                                 int &prop_depth,
                                 double &prop_value);

    virtual int compute_reverse_proposal_probability(int prop_idx,
                                                     int prop_depth,
                                                     double prop_value,
                                                     int &ii,
                                                     int &ij,
                                                     double &prop_parent_coeff,
                                                     double &prop_prob,
                                                     double &reverse_prob,
                                                     double &prior_prob);

    virtual int compute_likelihood(int prop_idx,
                                   double &proposed_likelihood,
                                   double &proposed_log_normalization);

    virtual int compute_acceptance(double proposed_likelihood,
                                   double proposed_log_normalization,
                                   double reverse_prob,
                                   double choose_prob,
                                   double prop_prob,
                                   double ratio,
                                   double prior_prob,
                                   bool &accept_proposal);

    virtual int communicate_acceptance(bool &accept_proposal);
}