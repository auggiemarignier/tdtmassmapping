#pragma once

#include "globalprop.hpp"

#include <mpi.h>

class Proposal
{
public:
    Proposal(GlobalProposal &global);
    Proposal(GlobalProposal &global, wavetree_perturb_t name);
    virtual ~Proposal()
    {
        delete[] propose_depth;
        delete[] accept_depth;
    };

    virtual int step();

    bool primary() const;

    std::string write_short_stats();

    std::string write_long_stats();

    void initialize_mpi(MPI_Comm communicator);

    GlobalProposal &global;
    wavetree_perturb_t name;
    int propose;
    int accept;

    int *propose_depth;
    int *accept_depth;

    MPI_Comm communicator;
    int mpi_size;
    int mpi_rank;

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
                                                   int &ij) = 0;

    virtual int communicate_proposal_location_and_value(int &prop_valid,
                                                        int &prop_idx,
                                                        int &prop_depth,
                                                        double &prop_value);

    virtual int propose_proposal(int &prop_valid,
                                 int &prop_idx,
                                 int &prop_depth,
                                 double &prop_value) = 0;

    virtual int compute_reverse_proposal_probability(int prop_idx,
                                                     int prop_depth,
                                                     double prop_value,
                                                     int &ii,
                                                     int &ij,
                                                     double &prop_parent_coeff,
                                                     double &prop_prob,
                                                     double &reverse_prob,
                                                     double &prior_prob);

    virtual int sub_reverse_proposal(int prop_idx,
                                     int prop_depth,
                                     double prop_value,
                                     int &ii,
                                     int &ij,
                                     double &prop_parent_coeff,
                                     double &prop_prob,
                                     double &reverse_prob,
                                     double &prior_prob) = 0;

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

    virtual double calculate_alpha(double proposed_likelihood,
                                   double proposed_log_normalization,
                                   double reverse_prob,
                                   double choose_prob,
                                   double prop_prob,
                                   double ratio,
                                   double prior_prob) = 0;

    virtual int communicate_acceptance(bool &accept_proposal);
};

class BirthProposal : public Proposal
{
public:
    BirthProposal(GlobalProposal &global)
        : Proposal(global, WT_PERTURB_BIRTH){};

private:
    int choose_proposal_location_and_value(int k,
                                           double &ratio,
                                           int &prop_depth,
                                           int &prop_idx,
                                           double &choose_prob,
                                           double &prop_value,
                                           double &prop_prob,
                                           int &prop_valid,
                                           double &prop_parent_coeff,
                                           int &ii,
                                           int &ij) override;

    int propose_proposal(int &prop_valid,
                         int &prop_idx,
                         int &prop_depth,
                         double &prop_value) override;

    int sub_reverse_proposal(int prop_idx,
                             int prop_depth,
                             double prop_value,
                             int &ii,
                             int &ij,
                             double &prop_parent_coeff,
                             double &prop_prob,
                             double &reverse_prob,
                             double &prior_prob) override;

    double calculate_alpha(double proposed_likelihood,
                           double proposed_log_normalization,
                           double reverse_prob,
                           double choose_prob,
                           double prop_prob,
                           double ratio,
                           double prior_prob) override;
};

class DeathProposal : public Proposal
{
public:
    DeathProposal(GlobalProposal &global)
        : Proposal(global, WT_PERTURB_DEATH){};

private:
    int choose_proposal_location_and_value(int k,
                                           double &ratio,
                                           int &prop_depth,
                                           int &prop_idx,
                                           double &choose_prob,
                                           double &prop_value,
                                           double &prop_prob,
                                           int &prop_valid,
                                           double &prop_parent_coeff,
                                           int &ii,
                                           int &ij) override;

    int choose_death_location(int k,
                              double &ratio,
                              int &death_depth,
                              int &death_idx,
                              double &choose_prob,
                              int &death_valid);

    int propose_proposal(int &prop_valid,
                         int &prop_idx,
                         int &prop_depth,
                         double &prop_value) override;

    int sub_reverse_proposal(int prop_idx,
                             int prop_depth,
                             double prop_value,
                             int &ii,
                             int &ij,
                             double &prop_parent_coeff,
                             double &prop_prob,
                             double &reverse_prob,
                             double &prior_prob) override;

    double calculate_alpha(double proposed_likelihood,
                           double proposed_log_normalization,
                           double reverse_prob,
                           double choose_prob,
                           double prop_prob,
                           double ratio,
                           double prior_prob) override;
};