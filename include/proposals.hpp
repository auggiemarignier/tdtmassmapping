#include <globalslice.hpp>
#include <birthslice.hpp>
#include <deathslice.hpp>
#include <valueslice.hpp>

#include <array>

using namespace std;

class GlobalSliceMM : public GlobalSlice
{
public:
    GlobalSliceMM(const char *filename);

    array<double, 256 * 256> inputdata;
    double mean;
    double var;
    double stddev;

private:
    void readdatafile(const char *filename);
    int n_obs;
};