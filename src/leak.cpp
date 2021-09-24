#include "mmobservations.hpp"
#include "rng.hpp"
#include "logging.hpp"

int main()
{
    uint imsize = 1024;
    Rng random(1);

    complexvector input;
    fftw_complex *inputfftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
    for (uint j = 0; j < imsize; j++)
    {
        input.emplace_back(random.normal(1.), random.normal(1.));
        inputfftw[j][0] = input[j].real();
        inputfftw[j][1] = input[j].imag();
    }
    mmobservations observations(32, 32);
    Logger::open_log(0);

    int i = 1;
    while (i < 1e6)
    {
        fftw_complex *outputfftw = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * imsize);
        observations.kaiser_squires(outputfftw, inputfftw);
        fftw_free(outputfftw);
        i++;
        if (i % 10000 == 0)
            INFO("%i iterations", i);
    }
}