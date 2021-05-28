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
    GlobalSliceMM(std::vector<double> _obs,
                  std::vector<double> _sigma,
                  const char *prior_file,
                  int degreex,
                  int degreey,
                  int seed,
                  int kmax,
                  int waveletxy);
    ~GlobalSliceMM();

    double likelihood(double &log_normalization);
    void accept();

    std::vector<double> inputdata;
    std::vector<double> stddev_v;
    mmobservations *observations;
    double mean;
    double var;
    double stddev;

private:
    void
    readdatafile(const char *filename);
    int n_obs;
};

class DeathSliceMM : public DeathSlice
{
public:
    DeathSliceMM(GlobalSliceMM &global);

    GlobalSliceMM &global;
};

class BirthSliceMM : public BirthSlice
{
public:
    BirthSliceMM(GlobalSliceMM &global);

    GlobalSliceMM &global;
};

class ValueSliceMM : public ValueSlice
{
public:
    ValueSliceMM(GlobalSliceMM &global);

    GlobalSliceMM &global;
};