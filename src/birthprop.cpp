extern "C"
{
#include "slog.h"
};
#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

#include "proposals.hpp"

int BirthProposal::choose_proposal_location_and_value(int k,
                                                      double &ratio,
                                                      int &birth_depth,
                                                      int &birth_idx,
                                                      double &choose_prob,
                                                      double &birth_value,
                                                      double &birth_prob,
                                                      int &birth_valid,
                                                      double &birth_parent_coeff,
                                                      int &ii,
                                                      int &ij)
{
    //
    // First determine coefficient to birth and new value
    //
    if (primary())
    {

        if (hnk_get_kplus1_ratio(global.hnk,
                                 global.treemaxdepth,
                                 k,
                                 &ratio) < 0)
        {
            ERROR("failed to get ratio for birth\n");
            return -1;
        }

        if (ratio <= 0.0)
        {
            ERROR("invalid ratio %d %d %f\n", global.treemaxdepth, k, ratio);
            return -1;
        }

        //
        // Note for birth we actually want the inverse of the ratio
        //
        ratio = 1.0 / ratio;

        if (wavetree2d_sub_choose_birth_global(global.wt,
                                               global.random.uniform(),
                                               global.treemaxdepth,
                                               &birth_depth,
                                               &birth_idx,
                                               &choose_prob) < 0)
        {
            /* Generally means full tree, might need to check this later. */
            ERROR("failed to choose birth\n");
            return true;
        }

        if (coefficient_histogram_propose_birth(global.coeff_hist, birth_idx) < 0)
        {
            ERROR("failed to update histogram for birth proposal\n");
            return -1;
        }

        if (wavetree2d_sub_2dindices(global.wt, birth_idx, &ii, &ij) < 0)
        {
            ERROR("failed to compute 2d indices for birth\n");
            return -1;
        }

        if (wavetree2d_sub_get_coeff(global.wt,
                                     wavetree2d_sub_parent_index(global.wt, birth_idx),
                                     &birth_parent_coeff) < 0)
        {
            ERROR("failed to get parent coefficient for birth (idx = %d)\n", birth_idx);
            return -1;
        }

        if (wavetree_pp_birth2d(global.proposal,
                                ii,
                                ij,
                                birth_depth,
                                global.treemaxdepth,
                                birth_parent_coeff,
                                &birth_value,
                                &birth_prob,
                                &birth_valid) < 0)
        {
            ERROR("BirthSlice::step: failed to do birth proposal\n");
            return -1;
        }
    }

    return 0;
}

int BirthProposal::propose_proposal(int &birth_valid,
                                    int &birth_idx,
                                    int &birth_depth,
                                    double &birth_value)
{
    //
    // BirthSlice the point (use birth from prior)
    //
    if (wavetree2d_sub_propose_birth(global.wt,
                                     birth_idx,
                                     birth_depth,
                                     birth_value) < 0)
    {
        ERROR("failed to birth point\n");
        return -1;
    }
    return 0;
}

int BirthProposal::compute_reverse_proposal_probability(int birth_idx,
                                                        int birth_depth,
                                                        double birth_value,
                                                        int &ii,
                                                        int &ij,
                                                        double &birth_parent_coeff,
                                                        double &birth_prob,
                                                        double &reverse_prob,
                                                        double &prior_prob)
{
    if (wavetree2d_sub_reverse_birth_global(global.wt,
                                            global.treemaxdepth,
                                            birth_depth,
                                            birth_idx,
                                            &reverse_prob) < 0)
    {
        ERROR("failed to reverse birth (global)\n");
        return -1;
    }
    return 0;
}