#include "globalprop.hpp"

extern "C"
{
#include "hnk_cartesian_nonsquare.h"
#include "slog.h"
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
      lambda(1.0),
      observations(_observations),
      model(nullptr),
      workspace(nullptr),
      mean_residual_n(0),
      residual(nullptr),
      mean_residual(nullptr),
      last_valid_residual(nullptr),
      residual_normed(nullptr),
      mean_residual_normed(nullptr),
      last_valid_residual_normed(nullptr),
      residuals_valid(false),
      residual_hist_bins(100),
      residual_hist_min(-5.0),
      residual_hist_max(5.0),
      residual_hist(nullptr),
      width(-1),
      height(-1),
      size(-1),
      ncoeff(-1),
      zoffset(nullptr),
      hierarchical(nullptr),
      current_likelihood(-1.0),
      coeff_hist(nullptr),
      random(seed),
      xywaveletf(nullptr),
      temperature(1.0),
      cov_n(-1)
{
    if (degreex < 0 || degreex >= 16 ||
        degreey < 0 || degreey >= 16)
    {
        throw WAVETOMO2DEXCEPTION("Degree(s) out of range: %d x %d x %d\n", degreex, degreey, degreez);
    }

    xywaveletf = wavelet_inverse_function_from_id(xywavelet);
    if (xywaveletf == nullptr)
    {
        throw WAVETOMO2DEXCEPTION("Invalid horizontal wavelet %d\n", xywavelet);
    }

    hierarchical = new independentgaussianhierarchicalmodel();
    hierarchical->setparameter(0, lambda);

    wt = wavetree2d_sub_create(degreex, degreey, 0.0);
    if (wt == NULL)
    {
        throw WAVETOMO2DEXCEPTION("Failed to create wavetree\n");
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

    residual_size = ntotal;
    residual = new double[residual_size];
    mean_residual = new double[residual_size];
    last_valid_residual = new double[residual_size];
    residual_normed = new double[residual_size];
    mean_residual_normed = new double[residual_size];
    last_valid_residual_normed = new double[residual_size];

    residual_hist = new int[residual_size * residual_hist_bins];

    reset_residuals();

    if (initial_model == NULL)
    {
        if (wavetree2d_sub_initialize(wt, 0.0) < 0)
        {
            throw WAVETOMO2DEXCEPTION("Failed to initialize wavetree\n");
        }
    }
    else
    {

        if (wavetree2d_sub_load(wt, initial_model) < 0)
        {
            throw WAVETOMO2DEXCEPTION("Failed to load initial model: %s\n", initial_model);
        }

        INFO("Loaded model with %d coefficients\n", wavetree2d_sub_coeff_count(wt));
    }

    //
    // Hnk Ratio
    //
    if (kmax > ncoeff)
    {
        INFO("Warning: kmax truncated to %d\n", ncoeff);
        kmax = ncoeff;
    }

    hnk = hnk_cartesian_nonsquare_2D_create_sub(degreex,
                                                degreey,
                                                kmax);
    if (hnk == NULL)
    {
        throw WAVETOMO2DEXCEPTION("Failed to create hnk table\n");
    }

    //
    // Chain History
    //
    ch = chain_history_create(CHAIN_STEPS);
    if (ch == nullptr)
    {
        throw WAVETOMO2DEXCEPTION("Failed to create chain history\n");
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
        throw WAVETOMO2DEXCEPTION("Failed to create coefficient histogram\n");
    }

    //
    // Create proposal structure.
    //
    if (prior_file != nullptr)
    {
        proposal = wavetree_pp_load(prior_file, seed, coeff_hist);
        if (proposal == NULL)
        {
            throw WAVETOMO2DEXCEPTION("Failed to load proposal file\n");
        }

        for (int i = 0; i < ncoeff; i++)
        {
            int coeffdepth = wavetree2d_sub_depthofindex(wt, i);

            int ii, ij;
            if (wavetree2d_sub_2dindices(wt, i, &ii, &ij) < 0)
            {
                throw WAVETOMO2DEXCEPTION("Failed to get 2d indices\n");
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
                throw WAVETOMO2DEXCEPTION("Failed to get coefficient range\n");
            }

            if (coefficient_histogram_set_range(coeff_hist,
                                                i,
                                                vmin,
                                                vmax) < 0)
            {
                throw WAVETOMO2DEXCEPTION("Failed to set coefficient histogram range\n");
            }
        }
    }
}

GlobalProposal::~GlobalProposal()
{
    delete hierarchical;
    delete observations;
    delete[] zoffset;

    delete[] model;
    delete[] workspace;

    delete[] residual;
    delete[] mean_residual;
    delete[] last_valid_residual;
    delete[] residual_normed;
    delete[] mean_residual_normed;
    delete[] last_valid_residual_normed;
    delete[] residual_hist;

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
    std::vector<double> image_model_v(image_model, image_model + size);
    return observations->single_frequency_likelihood(image_model_v,
                                                     hierarchical,
                                                     residual,
                                                     residual_normed,
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
        throw WAVETOMO2DEXCEPTION("Failed to map model to array\n");
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
        throw WAVETOMO2DEXCEPTION("Failed to do inverse transform on coefficients\n");
    }

    log_normalization = 0.0;
    std::vector<double> model_v(model, model + size);
    return observations->single_frequency_likelihood(model_v,
                                                     hierarchical,
                                                     residual,
                                                     residual_normed,
                                                     log_normalization);
}

void GlobalProposal::reset_residuals()
{
    mean_residual_n = 0;
    for (int i = 0; i < residual_size; i++)
    {
        residual[i] = 0.0;
        mean_residual[i] = 0.0;
        last_valid_residual[i] = 0.0;

        residual_normed[i] = 0.0;
        mean_residual_normed[i] = 0.0;
        last_valid_residual_normed[i] = 0.0;

        for (int j = 0; j < residual_hist_bins; j++)
        {
            residual_hist[i * residual_hist_bins + j] = 0;
        }
    }

    cov_n = 0;

    for (int i = 0; i < (int)cov_count.size(); i++)
    {
        int N = cov_count[i];

        for (int j = 0; j < N; j++)
        {
            cov_delta[i][j] = 0.0;
            cov_mu[i][j] = 0.0;
        }

        N = N * N;
        for (int j = 0; j < N; j++)
        {
            cov_sigma[i][j] = 0.0;
        }
    }
}

void GlobalProposal::invalidate_residuals()
{
    residuals_valid = false;
}

void GlobalProposal::accept()
{
    residuals_valid = true;
    for (int i = 0; i < residual_size; i++)
    {
        last_valid_residual[i] = residual[i];
        last_valid_residual_normed[i] = residual_normed[i];
    }

    update_residual_mean();
    update_residual_covariance();
}

void GlobalProposal::reject()
{
    update_residual_mean();
}

void GlobalProposal::update_residual_mean()
{
    mean_residual_n++;

    for (int i = 0; i < residual_size; i++)
    {

        double delta = last_valid_residual[i] - mean_residual[i];
        mean_residual[i] += delta / (double)(mean_residual_n);

        delta = last_valid_residual_normed[i] - mean_residual_normed[i];
        mean_residual_normed[i] += delta / (double)(mean_residual_n);

        int hi = (int)((last_valid_residual_normed[i] - residual_hist_min) / (residual_hist_max - residual_hist_min) * (double)residual_hist_bins);
        if (hi >= 0 && hi < residual_hist_bins)
        {
            residual_hist[i * residual_hist_bins + hi]++;
        }
    }
}

void GlobalProposal::update_residual_covariance()
{
    double *p = last_valid_residual;

    for (int k = 0; k < residual_size; k++)
    {

        cov_n++;

        for (int i = 0; i < (int)cov_count.size(); i++)
        {

            int N = cov_count[i];

            for (int j = 0; j < N; j++)
            {
                cov_delta[i][j] = (p[j] - cov_mu[i][j]) / (double)(cov_n);
                cov_mu[i][j] += cov_delta[i][j];
            }

            for (int j = 0; j < N; j++)
            {
                for (int l = j; l < N; l++)
                {

                    cov_sigma[i][j * N + l] +=
                        (double)(cov_n - 1) * cov_delta[i][j] * cov_delta[i][l] -
                        cov_sigma[i][j * N + l] / (double)(cov_n);
                }
            }

            p += N;
        }
    }
}

int GlobalProposal::get_residual_size() const
{
    return residual_size;
}

const double *
GlobalProposal::get_mean_residuals() const
{
    return mean_residual;
}

const double *
GlobalProposal::get_mean_normed_residuals() const
{
    return mean_residual_normed;
}

bool GlobalProposal::save_residuals(const char *filename)
{
    if (observations != nullptr)
    {
        return observations->save_residuals(filename, mean_residual, mean_residual_normed);
    }
    else
    {
        return true;
    }
}

bool GlobalProposal::save_residual_histogram(const char *filename) const
{
    if (observations != nullptr)
    {
        FILE *fp = fopen(filename, "w");
        if (fp == NULL)
        {
            ERROR("Failed to create file");
            return false;
        }

        fprintf(fp, "%d %d %f %f\n", residual_size, residual_hist_bins, residual_hist_min, residual_hist_max);
        for (int i = 0; i < residual_size; i++)
        {
            for (int j = 0; j < residual_hist_bins; j++)
            {
                fprintf(fp, "%d ", residual_hist[i * residual_hist_bins + j]);
            }
            fprintf(fp, "\n");
        }

        fclose(fp);
    }

    return true;
}

bool GlobalProposal::save_residual_covariance(const char *filename) const
{
    if (observations != nullptr)
    {
        FILE *fp = fopen(filename, "w");
        if (fp == NULL)
        {
            ERROR("Failed to create file\n");
            return false;
        }

        fprintf(fp, "%d\n", (int)cov_count.size());

        for (int i = 0; i < (int)cov_count.size(); i++)
        {

            int N = cov_count[i];

            fprintf(fp, "%d\n", N);
            for (int j = 0; j < N; j++)
            {
                fprintf(fp, "%.9g ", cov_mu[i][j]);
            }
            fprintf(fp, "\n");

            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    fprintf(fp, "%.9g ", cov_sigma[i][j * N + k]);
                }
                fprintf(fp, "\n");
            }
        }

        fclose(fp);
    }

    return true;
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
