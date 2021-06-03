extern "C"
{
#include "slog.h"
};
#include "wavetomo2dutil.hpp"
#include "wavetomo2dexception.hpp"

#include "proposals.hpp"

Proposal::Proposal(GlobalProposal &_global)
    : global(_global),
      name(NONE),
      propose(0),
      accept(0),
      propose_depth(new int[global.treemaxdepth + 1]),
      accept_depth(new int[global.treemaxdepth + 1]),
      communicator(MPI_COMM_NULL),
      mpi_size(-1),
      mpi_rank(-1)
{
    for (int i = 0; i <= global.treemaxdepth; i++)
    {
        propose_depth[i] = 0;
        accept_depth[i] = 0;
    }
}

Proposal::Proposal(GlobalProposal &_global, NAME _name)
    : global(_global),
      name(_name),
      propose(0),
      accept(0),
      propose_depth(new int[global.treemaxdepth + 1]),
      accept_depth(new int[global.treemaxdepth + 1]),
      communicator(MPI_COMM_NULL),
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

    int k = wavetree2d_sub_coeff_count(global.wt);

    if (wavetree2d_sub_set_invalid_perturbation(global.wt, name) < 0)
    {
        ERROR("failed to initialise birth perturbation\n");
        return -1;
    }

    if (true) // figure this out for the three proposals
    {
        double ratio;
        int prop_depth;
        int prop_idx;
        double choose_prob;
        double prop_value;
        int prop_valid = (name == 1 || name == 2) ? 1 : 0; // not sure why should be different for death proposal
        double prop_parent_coeff;
        double prop_prob;
        double reverse_prob;
        double prior_prob = 0.0;

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
                                               ii,
                                               ij) < 0)
        {
            return -1;
        }

        if (communicate_proposal_location_and_value(prop_valid,
                                                    prop_idx,
                                                    prop_depth,
                                                    prop_value) < 0)
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
                                   ratio,
                                   prior_prob,
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
                accept++;
                accept_depth[prop_depth]++;

                // save to histogram

                global.current_likelihood = proposed_likelihood;
                global.current_log_normalization = proposed_log_normalization;
                global.accept();

                return 1;
            }
            else
            {
                // save to histogram

                global.reject();

                return 0;
            }
        }
    }
}

std::string Proposal::write_short_stats()
{
    return mkformatstring("%s %6d/%6d %7.3f",
                          accept,
                          propose,
                          propose == 0 ? 0.0 : 100.0 * (double)accept / (double)propose);
}

std::string Proposal::write_long_stats()
{
    std::string s = mkformatstring("%s: %6d %7.3f:",
                                   propose,
                                   propose == 0 ? 0.0 : 100.0 * (double)accept / (double)propose);
    for (int i = 0; i <= global.treemaxdepth; i++)
    {
        s = s + mkformatstring("%7.3f ",
                               propose_depth[i] == 0 ? 0.0 : 100.0 * (double)accept_depth[i] / (double)propose_depth[i]);
    }

    return s;
}

void Proposal::initialize_mpi(MPI_Comm _communicator)
{
    MPI_Comm_dup(_communicator, &communicator);

    if (MPI_Comm_size(communicator, &mpi_size) != MPI_SUCCESS)
    {
        throw WAVETOMO2DEXCEPTION("MPI Failure\n");
    }
    if (MPI_Comm_rank(communicator, &mpi_rank) != MPI_SUCCESS)
    {
        throw WAVETOMO2DEXCEPTION("MPI Failure\n");
    }
}

bool Proposal::primary() const
{
    return (communicator == MPI_COMM_NULL || mpi_rank == 0);
}

int Proposal::communicate_acceptance(bool &accept_proposal)
{
    if (communicator != MPI_COMM_NULL)
    {
        int ta;
        if (mpi_rank == 0)
        {
            ta = (int)accept_proposal;
        }
        if (MPI_Bcast(&ta, 1, MPI_INT, 0, communicator) != MPI_SUCCESS)
        {
            throw WAVETOMO2DEXCEPTION("Failed to broadcast accept proposal\n");
        }
        if (mpi_rank != 0)
        {
            accept_proposal = (bool)ta;
        }
    }
    return 0;
}