extern "C"
{
#include "slog.h"
};
#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

#include "proposals.hpp"

int ValueProposal::step()
{
    propose++;

    if (wavetree2d_sub_set_invalid_perturbation(global.wt, WT_PERTURB_VALUE) < 0)
    {
        return -1;
    }

    int value_depth;
    int value_idx;
    double choose_prob;
    double value;
    int ii, ij;
    double proposed_likelihood;
    double proposed_log_normalization;
    double value_prior_ratio = 1.0; // This must be 1 before calling
    int valid_proposal = 0;
    int prior_errors = 0;

    if (choose_value_location_and_value(value_depth,
                                        value_idx,
                                        choose_prob,
                                        value,
                                        ii,
                                        ij,
                                        value_prior_ratio,
                                        prior_errors,
                                        valid_proposal) < 0)
    {
        return -1;
    }
    if (communicate_value_location_and_value(valid_proposal,
                                             value_idx,
                                             value_depth,
                                             value) < 0)
    {
        return -1;
    }
    if (valid_proposal)
    {
        propose_depth[value_depth]++;
        if (propose_value(valid_proposal,
                          value_idx,
                          value_depth,
                          value) < 0)
        {
            return -1;
        }
        if (compute_likelihood(value_idx, proposed_likelihood, proposed_log_normalization) < 0)
        {
            return -1;
        }
        bool accept_proposal = false;
        if (compute_acceptance(value_idx,
                               value_prior_ratio,
                               proposed_likelihood,
                               proposed_log_normalization,
                               accept_proposal) < 0)
        {
            return -1;
        }
        if (communicate_acceptance(accept_proposal) < 0)
        {
            return -1;
        }
        if (accept_proposal)
        {
            //
            // Accept
            //
            accept++;
            accept_depth[value_depth]++;

            if (coefficient_histogram_accept_value(global.coeff_hist, value_idx, value) < 0)
            {
                ERROR("failed to update histogram for value acceptance\n");
                return -1;
            }
            if (wavetree2d_sub_commit(global.wt) < 0)
            {
                ERROR("failed to commit value (%d %d)\n", valid_proposal, prior_errors);
                return -1;
            }

            global.current_likelihood = proposed_likelihood;
            global.current_log_normalization = proposed_log_normalization;
            global.accept();

            return 1;
        }
        else
        {
            //
            // Reject
            //
            if (coefficient_histogram_reject_value(global.coeff_hist, value_idx, value) < 0)
            {
                ERROR("failed to update histogram for value rejection\n");
                return -1;
            }
            if (wavetree2d_sub_undo(global.wt) < 0)
            {
                ERROR("failed to undo value\n");
                return -1;
            }

            global.reject();

            return 0;
        }
    }
    return 0;
}