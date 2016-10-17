/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Implements the non-linear transformations. Includes motion and
 * observation models to be used by the GaussCanonical class in constructing
 * joint Gaussian distributions.
 *************************************************************************/

#include <vector>
#include <iostream>
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "matops.hpp"
#include "vecset.hpp"
#include "v2vtransform.hpp"
#include "transforms.hpp"

using namespace std;

// MotionModel
std::vector< ColVector<double> > MotionModel::operator()(const ColVector<double>& x) const {
	
	ASSERT(x.size() == 6, "x has inconsistent dimensions, needs to be 6x1");
	std::vector< ColVector<double> > y(1); y[0].resize(6);
	
	y[0][0] = x[0] + x[1]*delta_t_;
	y[0][1] = x[1];
	y[0][2] = x[2] + x[3]*delta_t_;
	y[0][3] = x[3];
	y[0][4] = x[4] + x[5]*delta_t_ - (0.5)*(9.81)*pow(delta_t_, 2);
	y[0][5] = x[5] - (9.81)*(delta_t_);

	return y;
} // operator()

// SensorModel
std::vector< ColVector<double> > SensorModel::operator()(const ColVector<double>& x) const {
	
	ASSERT(x.size() == 6, "x has inconsistent dimensions, needs to be 6x1");
	std::vector< ColVector<double> > z(1); z[0].resize(2);

	z[0][0] = sqrt( pow(x[0], 2) + pow(x[2], 2) + pow(x[4], 2) );
	z[1][1] = sqrt( pow(x[1], 2) + pow(x[3], 2) + pow(x[5], 2) );

	return z;
} // operator()
