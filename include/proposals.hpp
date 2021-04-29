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

private:
    void readdatafile(const char *filename);
};