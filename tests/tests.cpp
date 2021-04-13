#include "gtest/gtest.h"

#include "proposals.hpp"
#include "mmobservations.hpp"

class GlobalTest : public ::testing::Test
{
protected:
    GlobalSliceMM *global;
    void SetUp() override
    {
        global = new GlobalSliceMM("../../../../data/Bolshoi_7_clean_256.txt");
    }

    void TearDown() override {
        delete global;
    }
};

TEST_F(GlobalTest, FileReadIn) {
    ASSERT_TRUE(global->file_read);
}

TEST(MMObsTest, MMPredsIdentity)
{
    const double model = 0.0;
    ASSERT_FALSE(model);

    const double model1 = 1.0;
    ASSERT_TRUE(model1);
}

int main(int argc, char **argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}