#include <wavetree_prior.h>
#include <wavetreepp.h>

/*
 * Coefficients generalised gaussian at each depth with variable beta
 */
wavetree_prior_t *
wavetree_prior_create_depth_full_generalised_gaussian(int ndepths,
                                                      double *va,
                                                      double *beta,
                                                      unsigned long int seed);

wavetree_pp_t *
load_wavetree_pp(const char *filename, unsigned long int seed, coefficient_histogram_t *histogram);

wavetree_prior_t *
load_prior_wavetree_pp(FILE *fp, unsigned long int seed);
