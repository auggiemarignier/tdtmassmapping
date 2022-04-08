#include <math.h>

#include "globalprop.hpp"
#include "logging.hpp"

extern "C"
{
#include "hnk_cartesian_nonsquare.h"
#include "wavetree_prior_dfggd.h"
};

static constexpr double LARGE_LIKELIHOOD = 1e99;
static constexpr int SUBTILE = 1;

const int CHAIN_STEPS = 1000000;

int GlobalProposal_coordtoindex(void *user, int i, int j, int k, int depth)
{
    wavetree2d_sub_t *wt = (wavetree2d_sub_t *)user;

    return wavetree2d_sub_from_2dindices(wt, i, j);
}

int GlobalProposal_indextocoord(void *user, int index, int *i, int *j, int *k, int *depth)
{
    wavetree2d_sub_t *wt = (wavetree2d_sub_t *)user;

    if (wavetree2d_sub_2dindices(wt, index, i, j) < 0)
    {
        return -1;
    }

    k = 0;
    *depth = wavetree2d_sub_depthofindex(wt, index);
    return 0;
}

GlobalProposal::GlobalProposal(Observations *_observations,
                               const char *initial_model,
                               const char *prior_file,
                               int _degreex,
                               int _degreey,
                               int seed,
                               int _kmax,
                               int xywavelet)
    : kmax(_kmax),
      treemaxdepth(-1),
      wt(nullptr),
      ch(nullptr),
      hnk(nullptr),
      proposal(nullptr),
      degreex(_degreex),
      degreey(_degreey),
      degreez(0),
      observations(_observations),
      model(nullptr),
      workspace(nullptr),
      width(-1),
      height(-1),
      size(-1),
      ncoeff(-1),
      zoffset(nullptr),
      current_likelihood(-1.0),
      current_prior(-1.0),
      current_unnormed_posterior(-1.0),
      coeff_hist(nullptr),
      random(seed),
      xywaveletf(nullptr),
      temperature(1.0),
      cov_n(-1)
{
    if (degreex < 0 || degreex >= 16 ||
        degreey < 0 || degreey >= 16)
    {
        throw ERROR("Degree(s) out of range: %d x %d x %d\n", degreex, degreey, degreez);
    }

    xywaveletf = wavelet_inverse_function_from_id(xywavelet);
    if (xywaveletf == nullptr)
    {
        throw ERROR("Invalid horizontal wavelet %d\n", xywavelet);
    }

    wt = wavetree2d_sub_create(degreex, degreey, 0.0);
    if (wt == NULL)
    {
        throw ERROR("Failed to create wavetree\n");
    }

    width = wavetree2d_sub_get_width(wt);
    height = wavetree2d_sub_get_height(wt);
    slicestride = width * height;
    size = wavetree2d_sub_get_size(wt);
    ncoeff = wavetree2d_sub_get_ncoeff(wt);
    treemaxdepth = wavetree2d_sub_maxdepth(wt);

    INFO("Image: %d x %d\n", width, height);

    model = new double[size];
    int workspacesize = width;
    if (height > workspacesize)
    {
        workspacesize = height;
    }
    workspace = new double[workspacesize];

    int ntotal = observations->n_obs;
    INFO("Data: %d total points\n", ntotal);

    if (initial_model == NULL)
    {
        if (wavetree2d_sub_initialize(wt, 0.0) < 0)
        {
            throw ERROR("Failed to initialize wavetree\n");
        }
    }
    else
    {

        if (wavetree2d_sub_load(wt, initial_model) < 0)
        {
            throw ERROR("Failed to load initial model: %s\n", initial_model);
        }

        INFO("Loaded model with %d coefficients\n", wavetree2d_sub_coeff_count(wt));
    }

    //
    // Hnk Ratio
    //
    if (kmax > ncoeff)
    {
        WARNING("kmax truncated to %d\n", ncoeff);
        kmax = ncoeff;
    }

    hnk = hnk_cartesian_nonsquare_2D_create_sub(degreex,
                                                degreey,
                                                kmax);
    if (hnk == NULL)
    {
        throw ERROR("Failed to create hnk table\n");
    }

    //
    // Chain History
    //
    ch = chain_history_create(CHAIN_STEPS);
    if (ch == nullptr)
    {
        throw ERROR("Failed to create chain history\n");
    }

    //
    // Initialse coeff histogram
    //
    coeff_hist = coefficient_histogram_create(ncoeff, 100, -1.0, 1.0,
                                              GlobalProposal_coordtoindex,
                                              GlobalProposal_indextocoord,
                                              wt);
    if (coeff_hist == NULL)
    {
        throw ERROR("Failed to create coefficient histogram\n");
    }

    //
    // Create proposal structure.
    //
    if (prior_file != nullptr)
    {
        proposal = load_wavetree_pp(prior_file, seed, coeff_hist);
        if (proposal == NULL)
        {
            throw ERROR("Failed to load proposal file\n");
        }

        for (int i = 0; i < ncoeff; i++)
        {
            int coeffdepth = wavetree2d_sub_depthofindex(wt, i);

            int ii, ij;
            if (wavetree2d_sub_2dindices(wt, i, &ii, &ij) < 0)
            {
                throw ERROR("Failed to get 2d indices\n");
            }

            double vmin, vmax;
            if (wavetree_pp_prior_range2d(proposal,
                                          ii,
                                          ij,
                                          coeffdepth,
                                          treemaxdepth,
                                          0.0,
                                          &vmin,
                                          &vmax) < 0)
            {
                throw ERROR("Failed to get coefficient range\n");
            }

            if (coefficient_histogram_set_range(coeff_hist,
                                                i,
                                                vmin,
                                                vmax) < 0)
            {
                throw ERROR("Failed to set coefficient histogram range\n");
            }
        }
    }
}

GlobalProposal::~GlobalProposal()
{
    delete[] zoffset;
    delete[] model;
    delete[] workspace;

    coefficient_histogram_destroy(coeff_hist);
    hnk_destroy(hnk);
    wavetree2d_sub_destroy(wt);
    wavetree_pp_destroy(proposal);
}

double
GlobalProposal::image_likelihood(const double *image_model,
                                 double &log_normalization)
{
    log_normalization = 0.0;
    complexvector image_model_v(image_model, image_model + size);
    return observations->single_frequency_likelihood(image_model_v,
                                                     log_normalization);
}

double
GlobalProposal::likelihood(double &log_normalization)
{
    //
    // Get tree model wavelet coefficients
    //
    memset(model, 0, sizeof(double) * size);
    if (wavetree2d_sub_map_to_array(wt, model, size) < 0)
    {
        throw ERROR("Failed to map model to array\n");
    }

    //
    // Inverse wavelet transform
    //
    if (generic_lift_inverse2d(model,
                               width,
                               height,
                               width,
                               workspace,
                               xywaveletf,
                               xywaveletf,
                               SUBTILE) < 0)
    {
        throw ERROR("Failed to do inverse transform on coefficients\n");
    }

    log_normalization = 0.0;
    complexvector model_v(model, model + size);
    return observations->single_frequency_likelihood(model_v,
                                                     log_normalization);
}

double GlobalProposal::prior()
{
    double p_k = 1. / kmax; // prior on k

    int k = wavetree2d_sub_coeff_count(wt);
    mpz_t a;
    mpz_init(a);
    hnk_get_hnk(hnk, treemaxdepth, k, a);
    uint n = mpz_get_ui(a);
    mpz_clear(a);
    double p_hnk = 1. / n; // prior on tree arrangement

    double log_p_x = wavetree2d_sub_logpriorprobability(wt, proposal); // prior on coefficient values

    return log(p_k) + log(p_hnk) + log_p_x;
}

double GlobalProposal::unnormed_posterior(const double &log_likelihood, const double &log_prior)
{
    return log_likelihood + log_prior;
}

void GlobalProposal::accept()
{
}

void GlobalProposal::reject()
{
}

generic_lift_inverse1d_step_t
GlobalProposal::wavelet_inverse_function_from_id(int id)
{
    switch (id)
    {
    case WAVELET_HAAR:
        return haar_lift_inverse1d_haar_step;

    case WAVELET_DAUB4:
        return daub4_dwt_inverse1d_daub4_step;

    case WAVELET_DAUB6:
        return daub6_dwt_inverse1d_daub6_step;

    case WAVELET_DAUB8:
        return daub8_dwt_inverse1d_daub8_step;

    case WAVELET_CDF97:
        return cdf97_lift_inverse1d_cdf97_step;

    case WAVELET_CDF97_PERIODIC:
        return cdf97_lift_periodic_inverse1d_cdf97_step;

    default:
        return nullptr;
    }
}

generic_lift_forward1d_step_t
GlobalProposal::wavelet_forward_function_from_id(int id)
{
    switch (id)
    {
    case WAVELET_HAAR:
        return haar_lift_forward1d_haar_step;

    case WAVELET_DAUB4:
        return daub4_dwt_forward1d_daub4_step;

    case WAVELET_DAUB6:
        return daub6_dwt_forward1d_daub6_step;

    case WAVELET_DAUB8:
        return daub8_dwt_forward1d_daub8_step;

    case WAVELET_CDF97:
        return cdf97_lift_forward1d_cdf97_step;

    case WAVELET_CDF97_PERIODIC:
        return cdf97_lift_periodic_forward1d_cdf97_step;

    default:
        return nullptr;
    }
}

void GlobalProposal::set_max_depth(int md)
{
    int maxd = wavetree2d_sub_maxdepth(wt);
    if (md > maxd)
    {
        md = maxd;
    }

    if (md < 1)
    {
        md = 1;
    }

    treemaxdepth = md;

    int dmaxk = hnk_get_maxk_at_h(hnk, treemaxdepth);

    printf("Depth maxk: %d\n", dmaxk);
    if (kmax > dmaxk)
    {
        kmax = dmaxk;
    }
    printf("New   maxk: %d\n", kmax);
}
