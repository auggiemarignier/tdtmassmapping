#include "mmobservations.hpp"

bool mmobservations::single_frequency_predictions(const double *model){
    if (!model) {
        INFO("!model");
        return false;
    }
    else
    {
        INFO("Returning same model");
        return true;
    }
}