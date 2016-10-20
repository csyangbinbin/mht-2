/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for the util.hpp's classes.
 *************************************************************************/

#include <iostream>
#include "gtest/gtest.h"
#include "genmat.hpp"
#include "emdw.hpp"
#include "gausscanonical.hpp"
#include "utils.hpp"

TEST (TestUtils, ReadMeanAndCov) {
	emdw::RVIds vars = {0, 1, 2, 3, 4, 5};
	ColVector<double> mu(6); mu *= 0;
	
	Matrix<double> S_mat = gLinear::zeros<double>(6, 6);
	for (unsigned i = 0; i < 6; i++) S_mat(i, i) = 0.25;

	rcptr<Factor> factor = uniqptr<Factor>(new GaussCanonical(vars, mu, S_mat));
	
	ColVector<double> mu_read = ReadMean(factor);
	Matrix<double> S_read = ReadCovariance(factor);

	EXPECT_EQ(mu, mu_read);
	EXPECT_EQ(S_mat, S_read);
}
