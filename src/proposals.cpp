#include <iostream>
#include <fstream>
#include "proposals.hpp"

#include "wavetomo2dexception.hpp"
extern "C"
{
#include "slog.h"
};

static constexpr int SUBTILE = 1;

GlobalSliceMM::GlobalSliceMM(
    const char *filename,
    const char *prior_file,
    int degreex,
    int degreey,
    int seed,
    int kmax,
    int waveletxy)
    : GlobalSlice(
          NULL,
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

    model = new double[size];
    int workspacesize = width;
    if (height > workspacesize)
    {
        workspacesize = height;
    }
    workspace = new double[workspacesize];
    residual_size = observations->n_obs;
    residual = new double[residual_size];
    mean_residual = new double[residual_size];
    last_valid_residual = new double[residual_size];
    residual_normed = new double[residual_size];
    mean_residual_normed = new double[residual_size];
    last_valid_residual_normed = new double[residual_size];
    residual_hist = new int[residual_size * residual_hist_bins];
}

GlobalSliceMM::GlobalSliceMM(
    std::vector<double> _obs,
    std::vector<double> _sigma,
    const char *prior_file,
    int degreex,
    int degreey,
    int seed,
    int kmax,
    int waveletxy)
    : GlobalSlice(
          NULL,
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
          true),
      inputdata(_obs),
      stddev_v(_sigma)
{
    observations = new mmobservations(inputdata, stddev_v);
    model = new double[size];
    int workspacesize = width;
    if (height > workspacesize)
    {
        workspacesize = height;
    }
    workspace = new double[workspacesize];
    residual_size = observations->n_obs;
    residual = new double[residual_size];
    mean_residual = new double[residual_size];
    last_valid_residual = new double[residual_size];
    residual_normed = new double[residual_size];
    mean_residual_normed = new double[residual_size];
    last_valid_residual_normed = new double[residual_size];
    residual_hist = new int[residual_size * residual_hist_bins];
}

GlobalSliceMM::~GlobalSliceMM()
{
    delete observations;

    delete[] model;
    delete[] workspace;

    delete[] residual;
    delete[] mean_residual;
    delete[] last_valid_residual;
    delete[] residual_normed;
    delete[] mean_residual_normed;
    delete[] last_valid_residual_normed;
    delete[] residual_hist;
}

double
GlobalSliceMM::likelihood(double &log_normalization)
{
    //
    // Get tree model wavelet coefficients
    //
    memset(model, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, model, size) < 0)
    {
        throw WAVETOMO2DEXCEPTION("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(
            model,
            width,
            height,
            width,
            workspace,
            xywaveletf,
            xywaveletf,
            SUBTILE) < 0)
    {
        throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    log_normalization = 0.0;
    std::vector<double> model_v(model, model + sizeof(double) * size);
    return observations->single_frequency_likelihood(
        model_v,
        hierarchical,
        residual,
        residual_normed,
        log_normalization);
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
        file.close();
    }
    else
    {
        throw WAVETOMO2DEXCEPTION("File not opened %s", filename);
    }
}