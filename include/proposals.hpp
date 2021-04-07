#include <globalslice.hpp>
#include <birthslice.hpp>
#include <deathslice.hpp>
#include <valueslice.hpp>

class GlobalSliceMM : public GlobalSlice
{
public:
    GlobalSliceMM();

private:
    int width;
};