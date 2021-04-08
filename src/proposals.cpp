#include <iostream>
#include <fstream>
#include "proposals.hpp"
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
    cout << "in GlobalSliceMM Constructor\n";
    cout << "width = " << width << "\n";
    cout << "Opening file " << filename << "\n";

    ifstream file (filename);
    string str;
    if (file.is_open())
    {
        for (int i = 0; i < 256 * 256; i++){
            file >> inputdata[i];
        }
    }
    else
    {
        cout << "File not opened";
    }
    file.close();
}
