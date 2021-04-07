#include <iostream>
#include "proposals.hpp"

GlobalSliceMM::GlobalSliceMM() : GlobalSlice("dummy",
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
    std::cout << "in GlobalSliceMM Constructor\n";
    std::cout << "width = " << width << "\n";
}
