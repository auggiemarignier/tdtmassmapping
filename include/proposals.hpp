#include <globalslice.hpp>
#include <birthslice.hpp>
#include <deathslice.hpp>
#include <valueslice.hpp>

class GlobalSliceMM : public GlobalSlice
{
public:
    GlobalSliceMM(const char *filename);

    double inputdata[256 * 256];

private:
    int width;
};