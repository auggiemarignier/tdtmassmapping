#include <iostream>
#include <fstream>
#include "proposals.hpp"

extern "C" {
  #include "slog.h"
};

using namespace std;

GlobalSliceMM::GlobalSliceMM(const char *filename) : GlobalSlice(NULL,
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
                                                                 true),
                                                     width(5)
{
    INFO("In GlobalSliceMM Constructor\n");
    INFO("width = %i \n", width);
    readdatafile(filename);
}

void GlobalSliceMM::readdatafile(const char *filename)
{
    INFO("Opening file %s \n", filename);

    ifstream file(filename);
    string str;
    if (file.is_open())
    {
        for (int i = 0; i < 256 * 256; i++)
        {
            file >> inputdata[i];
        }
        file_read = true;
    }
    else
    {
        INFO("File not opened");
    }
    file.close();
}