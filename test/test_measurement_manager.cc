/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for measurement_manager.hpp.
 *************************************************************************/
#include <iostream>
#include "gtest/gtest.h"
#include "anytype.hpp"
#include "emdw.hpp"
#include "measurement_manager.hpp"
#include "utils.hpp"

class MMTest : public testing::Test {

	protected:
		virtual void SetUp() {
			fileName_ = "data/test_case_6";
			N_ = 6;
		}

		virtual void TearDown() {}

	protected:
		std::string fileName_;
		unsigned N_;
};

TEST_F (MMTest, InitTest) {
	rcptr<MeasurementManager> mm = uniqptr<MeasurementManager>(new MeasurementManager(fileName_, N_));
	std::vector<ColVector<double>> points = mm->getSensorPoints(0, 50);
	std::vector<emdw::RVVals> vals = mm->getSensorMeasurements(0, 50); 
}
