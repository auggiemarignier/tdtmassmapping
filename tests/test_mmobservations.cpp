#include "gtest/gtest.h"
#include <vector>

#include "mmobservations.hpp"

class MMObsTest : public ::testing::Test
{
protected:
    void SetUp() override {
        std::vector<double> obs = {1, 2, 3, 4, 5};
        double sigma = 1.41;
        observations = new mmobservations(obs, sigma);
    }
    void TearDown() override {
        delete observations;
    }

    mmobservations *observations;
};

TEST_F(MMObsTest, MMPredsIdentity)
{
    std::vector<double> model = {1, 2, 3, 4, 5};
    ASSERT_EQ(model, observations->single_frequency_predictions(model));
}

TEST_F(MMObsTest, MMLikelihood)
{
    std::vector<double> model = {1, 2, 3, 4, 5};
    independentgaussianhierarchicalmodel *hmodel = nullptr;
    hmodel = new independentgaussianhierarchicalmodel();
    hmodel->setparameter(0, 1.0);
    double log_normalization = 0.0;
    double residual[5];
    double residual_norm[5];

    ASSERT_EQ(0, observations->single_frequency_likelihood(
                     model, hmodel, residual, residual_norm, log_normalization));

    std::vector<double> model2 = {0, 2, 3, 4, 5};
    ASSERT_NE(0, observations->single_frequency_likelihood(
                     model2, hmodel, residual, residual_norm, log_normalization));

}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}