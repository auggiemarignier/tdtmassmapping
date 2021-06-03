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