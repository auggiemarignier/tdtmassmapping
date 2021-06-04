#include "wavetree.h"
#include <string>

std::string enum_to_string(wavetree_perturb_t type)
{
    switch (type)
    {
    case WT_PERTURB_INVALID:
        return "Invalid";
    case WT_PERTURB_NONE:
        return "None";
    case WT_PERTURB_BIRTH:
        return "Birth";
    case WT_PERTURB_DEATH:
        return "Death";
    case WT_PERTURB_VALUE:
        return "Value";
    case WT_PERTURB_MOVE:
        return "Move";
    case WT_PERTURB_HIERARCHICAL:
        return "Hierarchical";
    case WT_PERTURB_PTEXCHANGE:
        return "PT Exchange";
    case WT_PERTURB_PTMODELEXCHANGE:
        return "PT Model Exchange";
    }
}
