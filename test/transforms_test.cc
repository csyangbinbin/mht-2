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
		const double kTimeStep = 0.5;
		uniqptr<V2VTransform> motion_model_;
		uniqptr<ColVector<double>> x_state_;
	
		virtual void SetUp() {
			motion_model_ = uniqptr<MotionModel>(new MotionModel(kTimeStep));
			x_state_ = uniqptr<ColVector<double>>(new ColVector<double>(6));
			for(unsigned i = 0; i < 6; i++) (*x_state_)[i] = 0;
		}

		virtual void TearDown() {}
}; // MotionModelTest

class SensorModelTest : public testing::Test {
	public:
		uniqptr<V2VTransform> sensor_model_;
	
		virtual void SetUp() {
			sensor_model_ = uniqptr<SensorModel>(new SensorModel());
		}

		virtual void TearDown() {}
};// SensorModelTest

// MotionModel Tests
TEST_F (MotionModelTest, IncorrectSize) {
	x_state_->resize(4);
	ASSERT_ANY_THROW((*motion_model_)(*x_state_));
} 

TEST_F (MotionModelTest, TranslateVector) {
	ColVector<double> expected_result(6);
	expected_result[0] = 0; expected_result[1] = 0;
	expected_result[2] = 0; expected_result[3] = 0;
	expected_result[4] = -1.22625; expected_result[5] = -4.905;
	
	std::vector<ColVector<double>> actual_result((*motion_model_)(*x_state_));
	EXPECT_EQ(expected_result, actual_result[0]);
}

// SensorModel Tests
TEST_F (SensorModelTest, PredictMeasurement) {
	std::vector<ColVector<double>> actual_result((*motion_model_)(*x_state_));
}
