#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gtest/gtest.h"

extern "C"
{
#include "wavetree_prior_dfggd.h"
}

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

    /*
   * Check generation works
   */
    for (int level = 0; level < 3; level++)
    {
        mean = 0.0;
        for (i = 0; i < 1000; i++)
        {
            ASSERT_GE(prior->sample(prior->user,
                                    1, 2, 3,
                                    level,
                                    3,
                                    0.0,
                                    &coeff),
                      0);

            mean = mean + coeff;

            ASSERT_GE(prior->prob(prior->user,
                                  1, 2, 3,
                                  level,
                                  3,
                                  0.0,
                                  coeff),
                      0.);
        }

        /*
        * Ensure that the mean is near 0
        */
        mean /= (double)1000;
        ASSERT_NEAR(fabs(mean), 0., 1e-1);
    }

    wavetree_prior_destroy(prior);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
