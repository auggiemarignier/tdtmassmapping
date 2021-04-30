#include <iostream>
#include <fstream>
#include "proposals.hpp"

extern "C"
{
#include "slog.h"
};

using namespace std;

GlobalSliceMM::GlobalSliceMM(
    const char *filename,
    const char *prior_file,
    int degreex,
    int degreey,
    int seed,
    int kmax,
    int waveletxy)
    : GlobalSlice(NULL,
                  NULL,
                  prior_file,
                  degreex,
                  degreey,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  seed,
                  kmax,
                  1.0,
                  true,
                  waveletxy,
                  true)
{
    readdatafile(filename);

    observations = new mmobservations(inputdata, stddev);
}

void GlobalSliceMM::readdatafile(const char *filename)
{
    INFO("Opening file %s \n", filename);

    ifstream file(filename);
    double element;
    if (file.is_open())
    {
        while (file >> element)
        {
            inputdata.push_back(element);
            mean += element;
            n_obs++;
        }
        mean /= n_obs;
        for (int i = 0; i < n_obs; i++)
        {
            var += (inputdata[i] - mean) * (inputdata[i] - mean);
        }
        var /= n_obs;
        stddev = sqrt(var);
    }
    else
    {
        INFO("File not opened");
    }
    file.close();
}