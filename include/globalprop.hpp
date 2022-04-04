#pragma once

#include <vector>
#include <string>
#include <set>
#include <cstring>

#include <gmp.h>

extern "C"
{
#include "hnk.h"
#include "wavetree2d_sub.h"
#include "wavetreepp.h"
#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"
#include "generic_lift.h"
};

#include "rng.hpp"
#include "mmobservations.hpp"

class GlobalProposal
{
public:
    enum
    {
        WAVELET_HAAR = 0,
        WAVELET_DAUB4 = 1,
        WAVELET_DAUB6 = 2,
        WAVELET_DAUB8 = 3,
        WAVELET_CDF97 = 4,
        WAVELET_CDF97_PERIODIC = 5,
        WAVELET_MAX = 5
    };

    GlobalProposal(Observations *observations,
                   const char *initial_model,
                   const char *prior_file,
                   int degreex,
                   int degreey,
                   int seed,
                   int kmax,
                   int waveletxy);
    ~GlobalProposal();

    double image_likelihood(const double *image_model,
                            double &log_normalization);
    double likelihood(double &log_normalization);
    double likelihood_mpi(double &log_normalization) { return 1.0; } //TODO: implement this
    double prior();
    void accept();
    void reject();
    void set_max_depth(int md);

    static generic_lift_inverse1d_step_t wavelet_inverse_function_from_id(int id);
    static generic_lift_forward1d_step_t wavelet_forward_function_from_id(int id);

    int kmax;
    int maxdepth;
    int treemaxdepth;

    bool linear;

    wavetree2d_sub_t *wt;
    chain_history_t *ch;

    hnk_t *hnk;
    wavetree_pp_t *proposal;

    int degreex;
    int degreey;
    int degreez;

    Observations *observations;

    double *model;
    double *workspace;

    int width;
    int height;
    int depth;
    int slicestride;
    int size;
    int ncoeff;
    double *zoffset;

    double current_likelihood;
    double current_log_normalization;
    double current_prior;

    coefficient_histogram_t *coeff_hist;

    Rng random;

    generic_lift_inverse1d_step_t xywaveletf;

    double temperature;

    int cov_n;
    std::vector<int> cov_count;
    std::vector<double *> cov_delta;
    std::vector<double *> cov_mu;
    std::vector<double *> cov_sigma;
};
