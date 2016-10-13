/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for the transform.h's classes.
 *************************************************************************/
#include <iostream>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "transforms.hpp"

class MotionModelTest : public testing::Test {
	public:
		rcptr<ColVector<double>> x_state;
		rcptr<V2VTransform> motion_model;
	
		virtual void SetUp() {
			x_state = uniqptr<ColVector<double>>(new ColVector<double>(6));
			motion_model = uniqptr<MotionModel> (new MotionModel());
		}

		virtual void TearDown() {}
};

TEST_F (MotionModelTest, TestSize) {
	EXPECT_EQ(6, x_state->size());
}
