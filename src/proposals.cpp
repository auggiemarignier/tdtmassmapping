extern "C"
{
#include "slog.h"
};

#include "proposals.hpp"

Proposal::Proposal(GlobalProposal &_global)
    : global(_global),
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

std::string Proposal::write_short_stats()
{
    return mkformatstring("Birth %6d/%6d %7.3f",
                          accept,
                          propose,
                          propose == 0 ? 0.0 : 100.0 * (double)accept / (double)propose);
}

std::string Proposal::write_long_stats()
{
    std::string s = mkformatstring("Birth: %6d %7.3f:",
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