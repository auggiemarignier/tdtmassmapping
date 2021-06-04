extern "C"
{
#include "slog.h"
};
#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

#include "proposals.hpp"

int DeathProposal::choose_proposal_location_and_value(int k,
                                                      double &ratio,
                                                      int &prop_depth,
                                                      int &prop_idx,
                                                      double &choose_prob,
                                                      double &prop_value,
                                                      double &prop_prob,
                                                      int &prop_valid,
                                                      double &prop_parent_coeff,
                                                      int &ii,
                                                      int &ij)
{
    return choose_death_location(k, ratio, prop_depth, prop_idx, choose_prob, prop_valid);
}

int DeathProposal::choose_death_location(int k,
                                         double &ratio,
                                         int &death_depth,
                                         int &death_idx,
                                         double &choose_prob,
                                         int &death_valid)
{
    if (primary())
    {
        //
        // Get tree structure ratio
        //

        if (hnk_get_kplus1_ratio(global.hnk,
                                 global.treemaxdepth,
                                 k - 1,
                                 &ratio) < 0)
        {
            ERROR("failed to get ratio for birth\n");
            return -1;
        }

        if (wavetree2d_sub_choose_death_global(global.wt,
                                               global.random.uniform(),
                                               global.treemaxdepth,
                                               &death_depth,
                                               &death_idx,
                                               &choose_prob) < 0)
        {
            death_valid = 0;
        }
    }
    return 0;
}

int DeathProposal::propose_proposal(int &death_valid,
                                    int &death_idx,
                                    int &death_depth,
                                    double &death_value)
{
    if (death_valid)
    {
        if (coefficient_histogram_propose_death(global.coeff_hist, death_idx) < 0)
        {
            ERROR("failed to update histogram for death proposal\n");
            return -1;
        }
        if (wavetree2d_sub_propose_death(global.wt,
                                         death_idx,
                                         death_depth,
                                         &death_value) < 0)
        {
            ERROR("failed to propose death\n");
            return -1;
        }
    }
    return 0;
}