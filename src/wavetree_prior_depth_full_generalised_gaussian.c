#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <slog.h>

#include "wavetree_prior_dfggd.h"

/*
 * Coefficients generalised gaussian at each depth
 */

struct depth_full_generalised_gaussian
{
    gsl_rng *rng;

    int ndepths;
    double *va;
    double *beta;
};

static int
range_depth_full_generalised_gaussian(void *user,
                                      int i,
                                      int j,
                                      int k,
                                      int level,
                                      int maxlevel,
                                      double parent_coeff,
                                      double *vmin,
                                      double *vmax)
{
    struct depth_full_generalised_gaussian *s;
    int d;

    s = (struct depth_full_generalised_gaussian *)user;

    d = level;
    if (d >= s->ndepths)
        d = s->ndepths - 1;

    /*
   * Just use 3 standard deviations for width
   */
    *vmax = s->va[d] * 3.0;
    *vmin = -(*vmax);

    return 0;
}

static int
sample_depth_full_generalised_gaussian(void *user,
                                       int i,
                                       int j,
                                       int k,
                                       int level,
                                       int maxlevel,
                                       double parent_coeff,
                                       double *coeff)
{
    struct depth_full_generalised_gaussian *s;
    int d;

    s = (struct depth_full_generalised_gaussian *)user;

    d = level;
    if (d >= s->ndepths)
        d = s->ndepths - 1;

    *coeff = gsl_ran_exppow(s->rng, s->va[d], s->beta[d]);

    return 0;
}

static double
prob_depth_full_generalised_gaussian(void *user,
                                     int i,
                                     int j,
                                     int k,
                                     int level,
                                     int maxlevel,
                                     double parent_coeff,
                                     double coeff)
{
    struct depth_full_generalised_gaussian *s;
    int d;
    double p;
    s = (struct depth_full_generalised_gaussian *)user;

    d = level;
    if (d >= s->ndepths)
        d = s->ndepths - 1;

    p = gsl_ran_exppow_pdf(coeff, s->va[d], s->beta[d]);
    return p;
}

static int
valid_depth_full_generalised_gaussian(void *user,
                                      int i,
                                      int j,
                                      int k,
                                      int level,
                                      int maxlevel,
                                      double parent_coeff,
                                      double coeff)
{
    /* All values are valid */
    return -1;
}

static int
setscale_depth_full_generalised_gaussian(void *user,
                                         double newscale,
                                         double *oldscale)
{
    return -1;
}

static int
destroy_depth_full_generalised_gaussian(void *user)
{
    struct depth_full_generalised_gaussian *s;

    s = (struct depth_full_generalised_gaussian *)user;

    gsl_rng_free(s->rng);
    free(s->va);
    free(s->beta);
    free(s);

    return 0;
}

wavetree_prior_t *
wavetree_prior_create_depth_full_generalised_gaussian(int ndepths,
                                                      double *va,
                                                      double *beta,
                                                      unsigned long int seed)
{
    wavetree_prior_t *w;
    struct depth_full_generalised_gaussian *s;
    int i;

    w = malloc(sizeof(wavetree_prior_t));
    if (w == NULL)
        return NULL;

    s = malloc(sizeof(struct depth_full_generalised_gaussian));
    if (s == NULL)
        return NULL;

    s->rng = gsl_rng_alloc(gsl_rng_taus);
    if (s->rng == NULL)
        return NULL;

    gsl_rng_set(s->rng, seed);

    s->ndepths = ndepths;

    s->va = malloc(sizeof(double) * ndepths);
    if (s->va == NULL)
        return NULL;

    s->beta = malloc(sizeof(double) * ndepths);
    if (s->beta == NULL)
        return NULL;

    for (i = 0; i < ndepths; i++)
    {
        s->va[i] = va[i];
        s->beta[i] = beta[i];
    }

    w->user = s;
    w->range = range_depth_full_generalised_gaussian;
    w->sample = sample_depth_full_generalised_gaussian;
    w->prob = prob_depth_full_generalised_gaussian;
    w->valid = valid_depth_full_generalised_gaussian;
    w->setscale = setscale_depth_full_generalised_gaussian;
    w->destroy = destroy_depth_full_generalised_gaussian;

    return w;
}

wavetree_pp_t *load_wavetree_pp(const char *filename, unsigned long int seed, coefficient_histogram_t *histogram)
{
    FILE *fp;

    wavetree_prior_t *prior;
    wavetree_bd_t *bd;
    wavetree_value_t *value;

    wavetree_pp_t *pp;

    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        return NULL;
    }

    prior = load_prior_wavetree_pp(fp, seed);
    if (prior == NULL)
    {
        return NULL;
    }

    bd = wavetree_pp_load_bd(fp, seed, prior);
    if (bd == NULL)
    {
        return NULL;
    }

    value = wavetree_pp_load_value(fp, seed, histogram);
    if (value == NULL)
    {
        return NULL;
    }

    pp = malloc(sizeof(wavetree_pp_t));
    if (pp == NULL)
    {
        return NULL;
    }

    pp->prior = prior;
    pp->bd = bd;
    pp->value = value;

    fclose(fp);

    return pp;
}

wavetree_prior_t *
load_prior_wavetree_pp(FILE *fp, unsigned long int seed)
{
    char buffer[256];
    double vmin[16];
    double vmax[16];
    int ndepths;
    double beta;
    int i;

    wavetree_prior_t *r;

    if (fp == NULL)
    {
        return NULL;
    }

    if (fscanf(fp, "%s\n", buffer) != 1)
    {
        ERROR("failed to read prior name");
        return NULL;
    }

    if (strcmp(buffer, "uniform") == 0)
    {

        if (fscanf(fp, "%lf %lf\n", &(vmin[0]), &(vmax[0])) != 2)
        {
            ERROR("failed to read uniform prior bounds");
            return NULL;
        }

        r = wavetree_prior_create_globally_uniform(vmin[0], vmax[0], seed);
        if (r == NULL)
        {
            ERROR("failed to create Globally Uniform prior");
            return NULL;
        }
    }
    else if (strcmp(buffer, "laplace") == 0)
    {

        if (fscanf(fp, "%lf\n", &beta) != 1)
        {
            ERROR("faield to read laplace prior width");
            return NULL;
        }

        r = wavetree_prior_create_globally_laplacian(beta, seed);
        if (r == NULL)
        {
            ERROR("failed to create Globally Laplacian prior");
            return NULL;
        }
    }
    else if (strcmp(buffer, "depthuniform") == 0)
    {

        if (fscanf(fp, "%d\n", &ndepths) != 1)
        {
            ERROR("failed to read n depths for Depth Uniform");
            return NULL;
        }

        if (ndepths < 0 || ndepths > 16)
        {
            ERROR("invalid ndepths %d", ndepths);
            return NULL;
        }

        for (i = 0; i < ndepths; i++)
        {
            if (fscanf(fp, "%lf %lf\n", &(vmin[i]), &(vmax[i])) != 2)
            {
                ERROR("failed to read std dev %d", i);
                return NULL;
            }
        }

        r = wavetree_prior_create_depth_uniform(ndepths,
                                                vmin,
                                                vmax,
                                                seed);
        if (r == NULL)
        {
            ERROR("failed to create Depth Uniform prior");
            return NULL;
        }

        /* } else if (strcmp(buffer, "depthlaplacian") == 0) { */

        /*   if (fscanf(fp, "%d\n", &ndepths) != 1) { */
        /*     ERROR("failed to read n depths for Depth Laplacian"); */
        /*     return NULL; */
        /*   } */

        /*   if (ndepths < 0 || ndepths > 16) { */
        /*     ERROR("invalid ndepths %d", ndepths); */
        /*     return NULL; */
        /*   } */

        /*   for (i = 0; i < ndepths; i ++) { */
        /*     if (fscanf(fp, "%lf\n", &(vmin[i])) != 1) { */
        /* 	ERROR("failed to read std dev %d", i); */
        /* 	return NULL; */
        /*     } */
        /*   } */

        /*   r = wavetree_prior_create_depth_laplacian(ndepths, vmin, seed); */

        /*   if (r == NULL) { */
        /*     ERROR("failed to create Depth Uniform prior"); */
        /*     return NULL; */
        /*   } */
    }
    else if (strcmp(buffer, "depthgeneralisedgaussian") == 0)
    {

        if (fscanf(fp, "%lf\n", &beta) != 1)
        {
            ERROR("failed to read prior beta");
            return NULL;
        }

        if (fscanf(fp, "%d\n", &ndepths) != 1)
        {
            ERROR("failed to read n depths for Generalise Gaussian");
            return NULL;
        }

        if (ndepths < 0 || ndepths > 16)
        {
            ERROR("invalid ndepths %d", ndepths);
            return NULL;
        }

        for (i = 0; i < ndepths; i++)
        {
            if (fscanf(fp, "%lf\n", &(vmin[i])) != 1)
            {
                ERROR("failed to read std dev %d", i);
                return NULL;
            }
        }

        r = wavetree_prior_create_depth_generalised_gaussian(ndepths, vmin, beta, seed);
        if (r == NULL)
        {
            ERROR("failed to create generalised gaussian");
            return NULL;
        }
    }
    else
    {
        ERROR("invalid bd proprosal name: %s", buffer);
        return NULL;
    }

    return r;
}