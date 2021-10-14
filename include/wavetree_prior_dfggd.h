#include <wavetree_prior.h>

/*
 * Coefficients generalised gaussian at each depth with variable beta
 */
wavetree_prior_t *
wavetree_prior_create_depth_full_generalised_gaussian(int ndepths,
                                                      double *va,
                                                      double *beta,
                                                      unsigned long int seed);