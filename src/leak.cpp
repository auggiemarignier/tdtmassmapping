#include "mmobservations.hpp"
#include "rng.hpp"
#include "logging.hpp"

int main()
{
    uint imsize = 1024;
    Rng random(1);

    complexvector input;
    for (uint j = 0; j < imsize; j++)
    {
        input.emplace_back(random.normal(1.), random.normal(1.));
    }
    mmobservations observations(32, 32);
    Logger::open_log(0);

    int i = 1;
    while (i < 1e6)
    {
        observations.single_frequency_predictions(input);
        i++;
        if (i % 10000 == 0)
            INFO("%i iterations", i);
    }
}