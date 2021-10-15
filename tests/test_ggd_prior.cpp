
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gtest/gtest.h"

#include "wavetree_prior_dfggd.h"

TEST(TestFGGD, Values)
{
    wavetree_prior_t *prior;
    int i;
    double coeff;
    double mean;

    int ndepths = 3;
    double va[3] = {1., 1., 1.};
    double beta[3] = {1., 1.5, 2.};

    prior = wavetree_prior_create_depth_full_generalised_gaussian(ndepths, va, beta, 1234);
    ASSERT_TRUE(prior != NULL);

#if 0
    /*
   * Check generation works
   */
    mean = 0.0;
    for (i = 0; i < 1000; i++)
    {
        ck_assert(prior->sample(prior->user,
                                1, 2, 3,
                                5,
                                6,
                                0.0,
                                &coeff) >= 0);

        ck_assert(coeff >= -1.0);
        ck_assert(coeff <= 1.0);

        mean = mean + coeff;

        ck_assert(prior->prob(prior->user,
                              1, 2, 3,
                              5,
                              6,
                              0.0,
                              coeff) == 0.5);

        ck_assert(prior->valid(prior->user,
                               1, 2, 3,
                               5,
                               6,
                               0.0,
                               coeff));
    }

    /*
   * Ensure that the mean is near 0
   */
    mean /= (double)1000;
    ck_assert(fabs(mean) < 0.1);

    /*
   * Check values outside of prior return invalid an 0 prob
   */
    ck_assert(prior->prob(prior->user,
                          1, 2, 3,
                          5,
                          6,
                          0.0,
                          2.0) == 0.0);
    ck_assert(prior->valid(prior->user,
                           1, 2, 3,
                           5,
                           6,
                           0.0,
                           2.0) == 0);

    ck_assert(prior->prob(prior->user,
                          1, 2, 3,
                          5,
                          6,
                          0.0,
                          -2.0) == 0.0);
    ck_assert(prior->valid(prior->user,
                           1, 2, 3,
                           5,
                           6,
                           0.0,
                           -2.0) == 0);
#endif

    wavetree_prior_destroy(prior);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
