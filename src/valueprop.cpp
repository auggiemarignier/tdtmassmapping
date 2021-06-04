extern "C"
{
#include "slog.h"
};
#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

#include "proposals.hpp"

int ValueProposal::choose_proposal_location_and_value(int k,
                                                      double &ratio,
                                                      int &prop_depth,
                                                      int &prop_idx,
                                                      double &choose_prob,
                                                      double &prop_value,
                                                      double &prop_prob,
                                                      int &prop_valid,
                                                      double &prop_parent_coeff,
                                                      double &value_prior_ratio,
                                                      int &prior_errors,
                                                      int &ii,
                                                      int &ij)
{
    return choose_value_location_and_value(prop_depth,
                                           prop_idx,
                                           choose_prob,
                                           prop_value,
                                           ii,
                                           ij,
                                           value_prior_ratio,
                                           prior_errors,
                                           prop_valid);
}

int ValueProposal::choose_value_location_and_value(int &value_depth,
                                                   int &value_idx,
                                                   double &choose_prob,
                                                   double &value,
                                                   int &ii,
                                                   int &ij,
                                                   double &value_prior_ratio,
                                                   int &prior_errors,
                                                   int &valid_proposal)
{
    if (primary())
    {

        if (wavetree2d_sub_choose_value_global(global.wt,
                                               global.random.uniform(),
                                               global.treemaxdepth,
                                               &value_depth,
                                               &value_idx,
                                               &choose_prob) < 0)
        {
            ERROR("failed to choose global value\n");
            return -1;
        }

        if (wavetree2d_sub_get_coeff(global.wt,
                                     value_idx,
                                     &value) < 0)
        {
            ERROR("failed to get coefficient value\n");
            return -1;
        }

        if (wavetree2d_sub_2dindices(global.wt, value_idx, &ii, &ij) < 0)
        {
            ERROR("failed to compute 2d indices for birth\n");
            return -1;
        }

        if (coefficient_histogram_propose_value(global.coeff_hist, value_idx) < 0)
        {
            ERROR("failed to update histogram for value proposal\n");
            return -1;
        }

        if (wavetree_pp_value_init(global.proposal) < 0)
        {
            ERROR("failed to initialize value proposal\n");
            return -1;
        }

        double value_parent_coeff = 0.0;

        if (wavetree_pp_propose_value2d(global.proposal,
                                        ii, ij,
                                        value_depth,
                                        global.maxdepth,
                                        value_parent_coeff,
                                        global.temperature,
                                        &value,
                                        &value_prior_ratio) < 0)
        {
            ERROR("failed to perturb value\n");
            return -1;
        }

        prior_errors = wavetree_pp_value_error_count(global.proposal);
        if (prior_errors < 0)
        {
            ERROR("failed to check errors\n");
            return -1;
        }

        if (prior_errors == 0)
        {
            valid_proposal = 1;
        }
    }

    return 0;
}

int ValueProposal::propose_proposal(int &valid_proposal,
                                    int &value_idx,
                                    int &value_depth,
                                    double &value)
{
    if (valid_proposal)
    {
        if (wavetree2d_sub_propose_value(global.wt, value_idx, value_depth, value) < 0)
        {
            ERROR("failed to propose value\n");
            return -1;
        }
    }
    return 0;
}

bool ValueProposal::k_valid(int &k) { return true; }
