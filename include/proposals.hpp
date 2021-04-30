#include <globalslice.hpp>
#include <birthslice.hpp>
#include <deathslice.hpp>
#include <valueslice.hpp>

#include <vector>

#include "mmobservations.hpp"

using namespace std;

class GlobalSliceMM : public GlobalSlice
{
public:
    GlobalSliceMM(const char *filename,
                  const char *prior_file,
                  int degreex,
                  int degreey,
                  int seed,
                  int kmax,
                  int waveletxy);

    std::vector<double> inputdata;
    mmobservations *observations;
    double mean;
    double var;
    double stddev;

private:
    void
    readdatafile(const char *filename);
    int n_obs;
};