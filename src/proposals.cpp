#include <math.h>

#include "logging.hpp"
#include "proposals.hpp"
#include "utils.hpp"

Proposal::Proposal(GlobalProposal &_global)
    : global(_global),
      name(WT_PERTURB_NONE),
      propose(0),
      accept(0),
      propose_depth(new int[global.treemaxdepth + 1]),
      accept_depth(new int[global.treemaxdepth + 1]),
      mpi_size(-1),
      mpi_rank(-1)
{
    for (int i = 0; i <= global.treemaxdepth; i++)
    {
        propose_depth[i] = 0;
        accept_depth[i] = 0;
    }
}

Proposal::Proposal(GlobalProposal &_global, wavetree_perturb_t _name)
    : global(_global),
      name(_name),
      propose(0),
      accept(0),
      propose_depth(new int[global.treemaxdepth + 1]),
      accept_depth(new int[global.treemaxdepth + 1]),
      mpi_size(-1),
      mpi_rank(-1)
{
    for (int i = 0; i <= global.treemaxdepth; i++)
    {
        propose_depth[i] = 0;
        accept_depth[i] = 0;
    }
}

int Proposal::step()
{
    propose++;

    if (wavetree2d_sub_set_invalid_perturbation(global.wt, name) < 0)
    {
        ERROR("failed to initialise %s perturbation\n", enum_to_string(name).c_str());
        return -1;
    }

    int k = wavetree2d_sub_coeff_count(global.wt);

    if (k_valid(k))
    {
        double ratio;
        int prop_depth;
        int prop_idx;
        double choose_prob;
        double prop_value;
        int prop_valid = (name == 1 || name == 3) ? 0 : 1; // not sure why should be different for death proposal
        double prop_parent_coeff;
        double prop_prob;
        double reverse_prob;
        double prior_prob = 0.0;
        int prior_errors = 0;

        double proposed_likelihood;
        double proposed_log_normalization;

        int ii;
        int ij;

        if (choose_proposal_location_and_value(k,
                                               ratio,
                                               prop_depth,
                                               prop_idx,
                                               choose_prob,
                                               prop_value,
                                               prop_prob,
                                               prop_valid,
                                               prop_parent_coeff,
                                               prior_errors,
                                               ii,
                                               ij) < 0)
        {
            return -1;
        }

        if (prop_valid)
        {
            propose_depth[prop_depth]++;

            if (propose_proposal(prop_valid,
                                 prop_idx,
                                 prop_depth,
                                 prop_value) < 0)
            {
                return -1;
            }

            if (compute_reverse_proposal_probability(prop_idx,
                                                     prop_depth,
                                                     prop_value,
                                                     ii,
                                                     ij,
                                                     prop_parent_coeff,
                                                     prop_prob,
                                                     reverse_prob,
                                                     prior_prob) < 0)
            {
                return -1;
            }

            if (compute_likelihood(prop_idx,
                                   proposed_likelihood,
                                   proposed_log_normalization) < 0)
            {
                return -1;
            }

            bool accept_proposal = false;

            if (compute_acceptance(proposed_likelihood,
                                   proposed_log_normalization,
                                   reverse_prob,
                                   choose_prob,
                                   prop_prob,
                                   prop_idx,
                                   ratio,
                                   prior_prob,
                                   accept_proposal) < 0)
            {
                return -1;
            }

            if (accept_proposal)
            {
                accept++;
                accept_depth[prop_depth]++;

                if (coefficient_histogram_accept(global.coeff_hist, prop_idx, prop_value) < 0)
                {
                    ERROR("failed to update histogram for %s acceptance\n", enum_to_string(name).c_str());
                    return -1;
                }

                if (coefficient_histogram_sample(global.coeff_hist, prop_idx, prop_value) < 0)
                {
                    ERROR("failed to update histogram");
                    return -1;
                }

                if (wavetree2d_sub_commit(global.wt) < 0)
                {
                    ERROR("failed to commit %s\n", enum_to_string(name).c_str());
                    return -1;
                }

                global.current_likelihood = proposed_likelihood;
                global.current_log_normalization = proposed_log_normalization;
                global.current_prior = global.prior();
                global.current_unnormed_posterior = global.unnormed_posterior(global.current_likelihood, global.current_prior);
                global.accept();

                return 1;
            }
            else
            {
                if (coefficient_histogram_reject(global.coeff_hist, prop_idx, prop_value) < 0)
                {
                    ERROR("failed to update histogram for %s rejection\n", enum_to_string(name).c_str());
                    return -1;
                }

                if (coefficient_histogram_sample(global.coeff_hist, prop_idx, prop_value) < 0)
                {
                    ERROR("failed to update histogram");
                    return -1;
                }

                if (wavetree2d_sub_undo(global.wt) < 0)
                {
                    ERROR("failed to undo %s\n", enum_to_string(name).c_str());
                    return -1;
                }

                global.reject();

                return 0;
            }
        }
        else
            return 0;
    }
    else
        return 0;
}

int Proposal::compute_likelihood(int prop_idx,
                                 double &proposed_likelihood,
                                 double &proposed_log_normalization)
{
    proposed_likelihood = global.likelihood(proposed_log_normalization);
    return 0;
}

std::string Proposal::write_short_stats()
{
    return mkformatstring("%s %6d/%6d %7.3f",
                          enum_to_string(name).c_str(),
                          accept,
                          propose,
                          propose == 0 ? 0.0 : 100.0 * (double)accept / (double)propose);
}

std::string Proposal::write_long_stats()
{
    std::string s = mkformatstring("%s: %6d %7.3f:",
                                   enum_to_string(name).c_str(),
                                   propose,
                                   propose == 0 ? 0.0 : 100.0 * (double)accept / (double)propose);
    for (int i = 0; i <= global.treemaxdepth; i++)
    {
        s = s + mkformatstring("%7.3f ",
                               propose_depth[i] == 0 ? 0.0 : 100.0 * (double)accept_depth[i] / (double)propose_depth[i]);
    }

    return s;
}

int Proposal::compute_reverse_proposal_probability(int prop_idx,
                                                   int prop_depth,
                                                   double prop_value,
                                                   int &ii,
                                                   int &ij,
                                                   double &prop_parent_coeff,
                                                   double &prop_prob,
                                                   double &reverse_prob,
                                                   double &prior_prob)
{
    if (sub_reverse_proposal(prop_idx,
                             prop_depth,
                             prop_value,
                             ii,
                             ij,
                             prop_parent_coeff,
                             prop_prob,
                             reverse_prob,
                             prior_prob) < 0)
    {
        ERROR("failed to reverse Proposal\n");
        return -1;
    }

    //
    // Compute the prior ratio
    //
    prior_prob = wavetree_pp_prior_probability2d(global.proposal,
                                                 ii,
                                                 ij,
                                                 prop_depth,
                                                 global.treemaxdepth,
                                                 prop_parent_coeff,
                                                 prop_value);
    return 0;
}

int Proposal::compute_acceptance(double proposed_likelihood,
                                 double proposed_log_normalization,
                                 double reverse_prob,
                                 double choose_prob,
                                 double prop_prob,
                                 int prop_idx,
                                 double ratio,
                                 double prior_prob,
                                 bool &accept_proposal)
{
    double u = log(global.random.uniform());
    double alpha = calculate_alpha(proposed_likelihood,
                                   proposed_log_normalization,
                                   reverse_prob,
                                   choose_prob,
                                   prop_prob,
                                   prop_idx,
                                   ratio,
                                   prior_prob);
    accept_proposal = u < alpha;
    return 0;
}