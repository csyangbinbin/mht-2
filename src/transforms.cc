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
	// Assert dimensional consistency
	ASSERT(x.size() == 6, "x has inconsistent dimensions, needs to be 6x1");
	std::vector< ColVector<double> > y(1); y[0].resize(6);
	
	// Shift target through the motion model
	y[0][0] = x[0] + x[1]*deltaT_;
	y[0][1] = x[1];
	y[0][2] = x[2] + x[3]*deltaT_;
	y[0][3] = x[3];
	y[0][4] = x[4] + x[5]*deltaT_ - (0.5)*(9.81)*pow(deltaT_, 2);
	y[0][5] = x[5] - (9.81)*(deltaT_);

	return y;
} // operator()

// SensorModel
std::vector< ColVector<double> > SensorModel::operator()(const ColVector<double>& x) const {
	// Assert dimensional consistency
	ASSERT(x.size() == 6, "x has inconsistent dimensions, needs to be 6x1");
	std::vector< ColVector<double> > z(1); z[0].resize(2); z[0] *= 0;
	
	// Determine the 'true' Range-Doppler measurements
	ColVector<double> relativePosition = x.slice(gLinear::gIndexRange(0, 4 , 2)) - sensorPosition_;
	z[0][0] = sqrt((relativePosition.transpose())*relativePosition);
	z[0][1] = (1.0/z[0][0])*( relativePosition.transpose()*(x.slice(gLinear::gIndexRange(1, 5, 2))) );

	// Adjust the range measurements
	z[0][0] += (z[0][1]*fc_*tp_)/(bw_);

	// Account for Doppler wrapping
	if (z[0][1] > vMax_) z[0][1] = -(z[0][1] - 2*vMax_);
	if (z[0][1] < -vMax_) z[0][1] = -(z[0][1] + 2*vMax_);

	return z;
} // operator()
