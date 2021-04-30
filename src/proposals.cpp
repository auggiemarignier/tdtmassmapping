#include <iostream>
#include <fstream>
#include "proposals.hpp"

extern "C"
{
#include "slog.h"
};

using namespace std;

GlobalSliceMM::GlobalSliceMM(
    const char *filename)
    : GlobalSlice(NULL,
                  NULL,
                  nullptr,
                  3,
                  3,
                  3,
                  0,
                  0,
                  -10.0,
                  10.0,
                  -10.0,
                  10.0,
                  1,
                  5,
                  1.0,
                  true,
                  4,
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