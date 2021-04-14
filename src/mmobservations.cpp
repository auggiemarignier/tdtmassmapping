#include <vector>

#include "mmobservations.hpp"
extern "C"
{
#include "slog.h"
};

std::vector<double> mmobservations::single_frequency_predictions(std::vector<double> model)
{
    std::vector<double> predictions = model;
    return predictions;
}