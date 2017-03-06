/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for canonical_gaussian_mixture.hpp.
 *************************************************************************/
#include <iostream>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "system_constants.hpp"
#include "transforms.hpp"
#include "utils.hpp"

class TransformsTest : public testing::Test {

	protected:

	protected:
		virtual void SetUp() {
		}

		virtual void TearDown() {
		}

	protected:
};

TEST_F (TransformsTest, MotionModelTest) {
	std::vector<rcptr<V2VTransform>> transform = initialiseMeasurementModels();
	ColVector<double> x(6);

	x[0] = -94.8463;
	x[1] = -12.8935;
	x[2] = 22.0664;
	x[3] = -0.0470;
	x[4] = 13.4131;
	x[5] = -8.4455;

	for (unsigned i = 0; i < 6; i++) std::cout << (transform[i]->operator()(x))[0] << std::endl;
}
