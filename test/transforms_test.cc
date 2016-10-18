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
#include "gausscanonical.hpp"
#include "transforms.hpp"

class MotionModelTest : public testing::Test {
	protected:
		virtual void SetUp() {
			motion_model_ = uniqptr<MotionModel>(new MotionModel(kTimeStep_));
			x_state_ = uniqptr<ColVector<double>>(new ColVector<double>(6));
			for(unsigned i = 0; i < 6; i++) (*x_state_)[i] = 0;
		}

		virtual void TearDown() {}

		const double kTimeStep_ = 0.5;
		uniqptr<V2VTransform> motion_model_;
		uniqptr<ColVector<double>> x_state_;
}; // MotionModelTest


class SensorModelTest : public testing::Test {
	protected:
		virtual void SetUp() {
			sensor_location_ = uniqptr<ColVector<double>>(new ColVector<double>(3));
			for(unsigned i = 0; i < 3; i++) (*sensor_location_)[i] = 0;
			
			x_state_ = uniqptr<ColVector<double>>(new ColVector<double>(6));
			for(unsigned i = 0; i < 6; i++) (*x_state_)[i] = 0;
			
			motion_model_ = uniqptr<MotionModel>(new MotionModel(kTimeStep_));
			sensor_model_ = uniqptr<SensorModel>(new SensorModel(*sensor_location_));
		}

		virtual void TearDown() {}

		const double kTimeStep_ = 0.5;
		uniqptr<ColVector<double>> sensor_location_;
		uniqptr<ColVector<double>> x_state_;
		uniqptr<V2VTransform> motion_model_;
		uniqptr<V2VTransform> sensor_model_;
};// SensorModelTest

// MotionModel Tests
TEST_F (MotionModelTest, IncorrectInputSize) {
	x_state_->resize(4);
	ASSERT_ANY_THROW((*motion_model_)(*x_state_));
} 

TEST_F (MotionModelTest, TranslateVector) {
	ColVector<double> expected_result(6);
	expected_result[0] = 0; expected_result[1] = 0;
	expected_result[2] = 0; expected_result[3] = 0;
	expected_result[4] = -1.22625; expected_result[5] = -4.905;
	
	std::vector<ColVector<double>> transform_result((*motion_model_)(*x_state_));
	EXPECT_EQ(expected_result, transform_result[0]);
}

TEST_F (MotionModelTest, CreateJointDistribution) {
	emdw::RVIds prior_vars = {0, 1, 2, 3, 4, 5};
	emdw::RVIds new_vars = {6, 7, 8, 9, 10, 11};
	
	ColVector<double> mu_0(6);
	mu_0[0] = 0; mu_0[1] = 72.86; mu_0[2] = 0;
	mu_0[3] = 0; mu_0[4] = 0; mu_0[5] = 14.04;

	Matrix<double> S_0 = gLinear::zeros<double>(6, 6);
	for (unsigned i = 0; i < 6; i++) S_0(i, i) = 1;

	Matrix<double> R_mat = gLinear::zeros<double>(6, 6);
	for (unsigned i = 0; i < 6; i++) R_mat(i, i) = 2.5;

	uniqptr<GaussCanonical> prior = uniqptr<GaussCanonical>(new GaussCanonical(prior_vars, S_0, mu_0, true));
	uniqptr<GaussCanonical> joint_distribution = 
		uniqptr<GaussCanonical>(new GaussCanonical(prior->copy(), *motion_model_, new_vars, R_mat, true));
	ColVector<double> mean = joint_distribution->getMean();
	
	ColVector<double> expected_mean(12); expected_mean *= 0;
	expected_mean[1] = 72.86; expected_mean[5] = 14.04;
	expected_mean[6] = 36.34; expected_mean[7]= 72.86;
	expected_mean[10] = 5.97375; expected_mean[11] = 9.135;

	for (unsigned i = 0; i < 12; i++) EXPECT_NEAR(expected_mean[i], mean[i], 0.2); 	
}

// SensorModel Tests
TEST_F (SensorModelTest, IncorrectInputSize) {
	x_state_->resize(4);
	ASSERT_ANY_THROW((*sensor_model_)(*x_state_));
}

TEST_F (SensorModelTest, PredictMeasurement) {
	ColVector<double> input_vector(6);
	input_vector[0] = 0; input_vector[1] = 0;
	input_vector[2] = 3; input_vector[3] = 3;
	input_vector[4] = 4; input_vector[5] = 4;

	ColVector<double> expected_result(2);
	expected_result[0] = 5; expected_result[1] = 5;
	
	std::vector<ColVector<double>> transform_result((*sensor_model_)(input_vector));
	EXPECT_EQ(expected_result, transform_result[0]);
}

TEST_F (SensorModelTest, CreateMeasurementDistribution) {
	emdw::RVIds prior_vars = {0, 1, 2, 3, 4, 5};
	emdw::RVIds markov_vars = {6, 7, 8, 9, 10, 11};
	emdw::RVIds measurement_vars = {12, 13};
	
	ColVector<double> mu_0(6);
	mu_0[0] = 0; mu_0[1] = 72.86; mu_0[2] = 0;
	mu_0[3] = 0; mu_0[4] = 0; mu_0[5] = 14.04;

	Matrix<double> S_0 = gLinear::zeros<double>(6, 6);
	for (unsigned i = 0; i < 6; i++) S_0(i, i) = 1;

	Matrix<double> R_mat = gLinear::zeros<double>(6, 6);
	for (unsigned i = 0; i < 6; i++) R_mat(i, i) = 2.5;

	Matrix<double> Q_mat = gLinear::zeros<double>(2, 2);
	Q_mat(0, 0) = 0.1;  Q_mat(1, 1) = 0.1;

	uniqptr<Factor> prior = uniqptr<GaussCanonical>(new GaussCanonical(prior_vars, S_0, mu_0, true));
	uniqptr<Factor> markov_distribution = 
		uniqptr<GaussCanonical>(new GaussCanonical(prior->copy(), *motion_model_, markov_vars, R_mat, true));
	uniqptr<GaussCanonical> measurement_distribution =
		uniqptr<GaussCanonical>(new GaussCanonical(markov_distribution->marginalize(markov_vars)->copy(), 
					*sensor_model_, measurement_vars, Q_mat, true));
	
	ColVector<double> mean = measurement_distribution->getMean();
	std::cout << mean << std::endl;
}
