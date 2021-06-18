#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <getopt.h>

extern "C"
{
#include "chain_history.h"

#include "cdf97_lift.h"
#include "cdf97_lift_periodic.h"
#include "haar_lift.h"
#include "daub4_dwt.h"
#include "daub6_dwt.h"
#include "daub8_dwt.h"
#include "generic_lift.h"
};

#include "globalprop.cpp"

struct user_data
{
    int thincounter;
    int thin;
    int skip;

    int counter;

    int degree_max;
    int width;
    int height;
    int size;

    double *mean;
    double *variance;

    int **hist;
    int bins;
    double vmin;
    double vmax;

    double *model;
    double *workspace;

    generic_lift_inverse1d_step_t hwaveletf;

    wavetree2d_sub_t *wt;

    FILE *fp_out;
};

static const double CREDIBLE_INTERVAL = 0.95;

static int process(int i,
                   void *user,
                   const chain_history_change_t *step,
                   const multiset_int_double_t *S_v);

static int histogram_index(double v, double vmin, double vmax, int bins);

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins);
static double median_from_histogram(int *hist, double vmin, double vmax, int bins);
static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);
static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop);
static double hpd_from_histogram(int *hist, double vmin, double vmax, int bins, double hpd_interval, double &hpd_min, double &hpd_max);

int main(int argc, char *argv[])
{
    int c;
    int option_index;

    chain_history_t *ch;

    const char *input_file;
    const char *output_file;
    const char *stddev_file;
    const char *likelihood_file;

    int degree_x;
    int degree_y;

    int thin;
    int skip;
    int maxsteps;

    const char *mode_file;
    const char *median_file;
    const char *credible_min;
    const char *credible_max;
    const char *histogram;
    const char *hpd_min;
    const char *hpd_max;
    const char *hpd_range;

    int bins;
    double vmin;
    double vmax;

    FILE *fp_in;
    FILE *fp_out;

    struct user_data data;
    multiset_int_double_t *S_v;

    int i;
    int j;

    int credible_drop;

    int waveleth;

    /*
   * Default values
   */
    fp_in = NULL;
    fp_out = NULL;
    degree_x = 8;
    degree_y = 8;

    input_file = "../outputs/ch.dat";
    output_file = "../outputs/mean.txt";
    stddev_file = NULL;
    likelihood_file = "../outputs/likelihood.txt";

    mode_file = NULL;
    median_file = NULL;
    credible_min = NULL;
    credible_max = NULL;
    histogram = NULL;

    hpd_min = NULL;
    hpd_max = NULL;
    hpd_range = NULL;

    bins = 1000;
    vmin = 2.0;
    vmax = 4.0;

    thin = 0;
    skip = 0;

    maxsteps = 1000000;

    waveleth = 4;

    vmin = 0.001;
    vmax = 1.0;

    ch = chain_history_create(maxsteps);
    if (ch == NULL)
    {
        fprintf(stderr, "error: failed to create chain history\n");
        return -1;
    }

    data.thincounter = 0;
    data.thin = thin;
    data.skip = skip;

    data.counter = 0;

    data.wt = wavetree2d_sub_create(degree_x, degree_y, 0.0);
    if (data.wt == nullptr)
    {
        fprintf(stderr, "error: failed to create wavetree\n");
        return -1;
    }

    data.width = wavetree2d_sub_get_width(data.wt);
    data.height = wavetree2d_sub_get_height(data.wt);
    data.size = wavetree2d_sub_get_size(data.wt);

    printf("Image: %d x %d\n", data.width, data.height);

    data.mean = new double[data.size];
    memset(data.mean, 0, sizeof(double) * data.size);

    data.bins = bins;
    data.vmin = vmin;
    data.vmax = vmax;
    data.hist = new int *[data.size];
    for (i = 0; i < data.size; i++)
    {
        data.hist[i] = new int[data.bins];
        memset(data.hist[i], 0, sizeof(int) * data.bins);
    }

    data.variance = new double[data.size];
    memset(data.variance, 0, sizeof(double) * data.size);

    data.model = new double[data.size];

    int workspacesize = data.width;
    if (data.height > data.width)
    {
        workspacesize = data.height;
    }
    data.workspace = new double[workspacesize];

    data.hwaveletf = GlobalProposal::wavelet_inverse_function_from_id(waveleth);

    S_v = multiset_int_double_create();
    if (S_v == NULL)
    {
        fprintf(stderr, "error: failed to create multiset\n");
        return -1;
    }

    fp_in = fopen(input_file, "r");
    if (fp_in == NULL)
    {
        fprintf(stderr, "error: failed to open input file\n");
        return -1;
    }

    //
    // Likelihood output
    //
    fp_out = fopen(likelihood_file, "w");
    if (fp_out == NULL)
    {
        fprintf(stderr, "error: failed to open likelihood file\n");
        return -1;
    }
    data.fp_out = fp_out;

    /*
   * Process the chain history
   */
    while (!feof(fp_in))
    {

        if (chain_history_read(ch,
                               (ch_read_t)fread,
                               fp_in) < 0)
        {
            if (feof(fp_in))
            {
                break;
            }

            fprintf(stderr, "error: failed to read chain history\n");
            return -1;
        }

        if (chain_history_replay(ch,
                                 S_v,
                                 (chain_history_replay_function_t)process,
                                 &data) < 0)
        {
            fprintf(stderr, "error: failed to replay\n");
            return -1;
        }
    }
    printf("%d records\n", data.counter);
    printf("%d total\n", data.thincounter);
    fclose(fp_in);
    fclose(fp_out);

    /*
   * Mean output
   */
    fp_out = fopen(output_file, "w");
    if (fp_out == NULL)
    {
        fprintf(stderr, "error: failed to open mean file\n");
        return -1;
    }

    for (j = 0; j < data.height; j++)
    {
        for (i = 0; i < data.width; i++)
        {
            fprintf(fp_out, "%10.6f ", data.mean[j * data.width + i]);
        }
        fprintf(fp_out, "\n");
    }

    fclose(fp_out);

    for (i = 0; i < data.size; i++)
    {
        data.variance[i] /= (double)(data.counter - 1);
    }

    /*
   * Std. Deviation output
   */
    if (stddev_file != NULL)
    {
        fp_out = fopen(stddev_file, "w");
        if (fp_out == NULL)
        {
            fprintf(stderr, "error: failed to open stddev file\n");
            return -1;
        }

        for (j = 0; j < data.height; j++)
        {
            for (i = 0; i < data.width; i++)
            {
                fprintf(fp_out, "%10.6f ", sqrt(data.variance[j * data.width + i]));
            }
            fprintf(fp_out, "\n");
        }

        fclose(fp_out);
    }

    /*
   * Mode output
   */
    if (mode_file != NULL)
    {
        fp_out = fopen(mode_file, "w");
        if (fp_out == NULL)
        {
            fprintf(stderr, "error: failed to open mode file\n");
            return -1;
        }

        for (j = 0; j < data.height; j++)
        {
            for (i = 0; i < data.width; i++)
            {
                fprintf(fp_out, "%10.6f ", mode_from_histogram(data.hist[j * data.width + i], data.vmin, data.vmax, data.bins));
            }
            fprintf(fp_out, "\n");
        }

        fclose(fp_out);
    }

    /*
   * Median output
   */
    if (median_file != NULL)
    {
        fp_out = fopen(median_file, "w");
        if (fp_out == NULL)
        {
            fprintf(stderr, "error: failed to open median file\n");
            return -1;
        }

        for (j = 0; j < data.height; j++)
        {
            for (i = 0; i < data.width; i++)
            {
                fprintf(fp_out, "%10.6f ", median_from_histogram(data.hist[j * data.width + i], data.vmin, data.vmax, data.bins));
            }
            fprintf(fp_out, "\n");
        }

        fclose(fp_out);
    }

    /*
   * Credible Min
   */
    credible_drop = (int)(((double)data.counter * (1.0 - CREDIBLE_INTERVAL)) / 2.0);

    if (credible_min != NULL)
    {
        fp_out = fopen(credible_min, "w");
        if (fp_out == NULL)
        {
            fprintf(stderr, "error: failed to open credible min file\n");
            return -1;
        }

        for (j = 0; j < data.height; j++)
        {
            for (i = 0; i < data.width; i++)
            {
                fprintf(fp_out, "%10.6f ", head_from_histogram(data.hist[j * data.width + i], data.vmin, data.vmax, data.bins, credible_drop));
            }
            fprintf(fp_out, "\n");
        }

        fclose(fp_out);
    }

    /*
   * Credible Max
   */
    if (credible_max != NULL)
    {
        fp_out = fopen(credible_max, "w");
        if (fp_out == NULL)
        {
            fprintf(stderr, "error: failed to open credible max file\n");
            return -1;
        }

        for (j = 0; j < data.height; j++)
        {
            for (i = 0; i < data.width; i++)
            {
                fprintf(fp_out, "%10.6f ", tail_from_histogram(data.hist[j * data.width + i], data.vmin, data.vmax, data.bins, credible_drop));
            }
            fprintf(fp_out, "\n");
        }

        fclose(fp_out);
    }

    if (histogram != NULL)
    {
        fp_out = fopen(histogram, "w");
        if (fp_out == NULL)
        {
            fprintf(stderr, "error: failed to open histogram file\n");
            return -1;
        }

        fprintf(fp_out, "%d %d\n", data.size, data.bins);
        fprintf(fp_out, "%.6f %.6f\n", data.vmin, data.vmax);

        for (j = 0; j < data.size; j++)
        {
            for (i = 0; i < data.bins; i++)
            {

                fprintf(fp_out, "%d ", data.hist[j][i]);
            }
            fprintf(fp_out, "\n");
        }

        fclose(fp_out);
    }

    if (hpd_range != NULL || hpd_min != NULL || hpd_max != NULL)
    {

        FILE *fpr = NULL;
        FILE *fpmin = NULL;
        FILE *fpmax = NULL;

        if (hpd_range != NULL)
        {
            fpr = fopen(hpd_range, "w");
            if (fpr == NULL)
            {
                fprintf(stderr, "error: failed to open hpd range file\n");
                return -1;
            }
        }

        if (hpd_min != NULL)
        {
            fpmin = fopen(hpd_min, "w");
            if (fpmin == NULL)
            {
                fprintf(stderr, "error: failed to open hpd min file\n");
                return -1;
            }
        }

        if (hpd_max != NULL)
        {
            fpmax = fopen(hpd_max, "w");
            if (fpmax == NULL)
            {
                fprintf(stderr, "error: failed to open hpd max file\n");
                return -1;
            }
        }

        for (j = 0; j < data.height; j++)
        {
            for (i = 0; i < data.width; i++)
            {
                double hmin, hmax;

                double hrange = hpd_from_histogram(data.hist[j * data.width + i],
                                                   data.vmin, data.vmax, data.bins,
                                                   CREDIBLE_INTERVAL,
                                                   hmin, hmax);

                if (fpr != NULL)
                {
                    fprintf(fpr, "%10.6f ", hrange);
                }

                if (fpmin != NULL)
                {
                    fprintf(fpmin, "%10.6f ", hmin);
                }

                if (fpmax != NULL)
                {
                    fprintf(fpmax, "%10.6f ", hmax);
                }
            }

            if (fpr != NULL)
            {
                fprintf(fpr, "\n");
            }

            if (fpmin != NULL)
            {
                fprintf(fpmin, "\n");
            }

            if (fpmax != NULL)
            {
                fprintf(fpmax, "\n");
            }
        }

        if (fpr != NULL)
        {
            fclose(fpr);
        }

        if (fpmin != NULL)
        {
            fclose(fpmin);
        }

        if (fpmax != NULL)
        {
            fclose(fpmax);
        }
    }

    chain_history_destroy(ch);
    multiset_int_double_destroy(S_v);

    delete[] data.mean;
    delete[] data.variance;
    delete[] data.model;
    delete[] data.workspace;

    return 0;
}

static int process(int stepi,
                   void *user,
                   const chain_history_change_t *step,
                   const multiset_int_double_t *S_v)
{
    struct user_data *d = (struct user_data *)user;
    double delta;
    int i;
    int hi;

    if ((d->thincounter >= d->skip) && (d->thin <= 1 || ((d->thincounter - d->skip) % d->thin) == 0))
    {
        fprintf(d->fp_out, "%.6f\n", step->header.likelihood);

        memset(d->model, 0, sizeof(double) * d->size);

        if (wavetree2d_sub_set_from_S_v(d->wt, S_v) < 0)
        {
            fprintf(stderr, "process: failed to set wavetree (sub)\n");
            return -1;
        }

        if (wavetree2d_sub_map_to_array(d->wt, d->model, d->size) < 0)
        {
            fprintf(stderr, "process: failed to map to array\n");
            return -1;
        }

        if (generic_lift_inverse2d(d->model,
                                   d->width,
                                   d->height,
                                   d->width,
                                   d->workspace,
                                   d->hwaveletf,
                                   d->hwaveletf,
                                   SUBTILE) < 0)
        {
            throw ERROR("Failed to do inverse transform on coefficients\n");
        }

        /*
     * Update mean/variance calculation
     */
        d->counter++;
        for (i = 0; i < d->size; i++)
        {

            delta = d->model[i] - d->mean[i];

            d->mean[i] += delta / (double)(d->counter);
            d->variance[i] += delta * (d->model[i] - d->mean[i]);
        }

        /*
     * Update the histogram
     */
        for (i = 0; i < d->size; i++)
        {
            hi = histogram_index(d->model[i], d->vmin, d->vmax, d->bins);

            d->hist[i][hi]++;
        }
    }
    d->thincounter++;

    return 0;
}

static int histogram_index(double v, double vmin, double vmax, int bins)
{
    int i;

    i = (int)((double)bins * (v - vmin) / (vmax - vmin));

    if (i < 0)
    {
        return 0;
    }

    if (i > (bins - 1))
    {
        return bins - 1;
    }

    return i;
}

static double mode_from_histogram(int *hist, double vmin, double vmax, int bins)
{
    int i;
    int m;
    int mi;

    m = 0;
    mi = -1;

    for (i = 0; i < bins; i++)
    {
        if (hist[i] > m)
        {
            m = hist[i];
            mi = i;
        }
    }

    if (mi < 0)
    {
        return 0.0;
    }

    return ((double)mi + 0.5) / (double)bins * (vmax - vmin) + vmin;
}

static double median_from_histogram(int *hist, double vmin, double vmax, int bins)
{
    int i;
    int j;
    int ci;
    int cj;

    i = 0;
    j = bins - 1;
    ci = 0;
    cj = 0;

    while (i != j)
    {
        if (ci < cj)
        {
            ci += hist[i];
            i++;
        }
        else
        {
            cj += hist[j];
            j--;
        }
    }

    return ((double)i + 0.5) / (double)bins * (vmax - vmin) + vmin;
}

static double head_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
    int i;
    int ci;

    i = 0;
    ci = 0;
    while (i < bins && ci < drop)
    {
        if (hist[i] + ci >= drop)
        {
            break;
        }

        ci += hist[i];
        i++;
    }

    return ((double)i + 0.5) / (double)bins * (vmax - vmin) + vmin;
}

static double tail_from_histogram(int *hist, double vmin, double vmax, int bins, int drop)
{
    int i;
    int ci;

    i = bins - 1;
    ci = 0;
    while (i > 0 && ci < drop)
    {
        if (hist[i] + ci >= drop)
        {
            break;
        }

        ci += hist[i];
        i--;
    }

    return ((double)i + 0.5) / (double)bins * (vmax - vmin) + vmin;
}

static double hpd_from_histogram(int *hist, double vmin, double vmax, int bins, double hpd_interval, double &hpd_min, double &hpd_max)
{
    int sum;
    int mincount;

    //
    // First count number of samples
    //

    sum = 0;
    for (int i = 0; i < bins; i++)
    {
        sum += hist[i];
    }

    mincount = (int)(hpd_interval * (double)sum);

    //
    // Now brute force search for minimum hpd with edges at
    //
    double minwidth = vmax - vmin;
    double minleft = vmin;
    double minright = vmax;

    for (int i = 0; i < bins; i++)
    {

        double left = vmin + (double)i / (double)bins * (vmax - vmin);
        int j = i + 1;

        int count = hist[i];
        while (j < bins && count < mincount)
        {
            count += hist[j];
            j++;
        }

        if (count >= mincount)
        {

            double right = vmin + (double)j / (double)bins * (vmax - vmin);
            if (right - left < minwidth)
            {

                minwidth = right - left;
                minleft = left;
                minright = right;
            }
        }
    }

    hpd_min = minleft;
    hpd_max = minright;

    return minwidth;
}
