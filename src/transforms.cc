/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Implements the non-linear transformations. Includes motion and
 * observation models to be used by the GaussConical class in constructing
 * joint Gaussian distributions.
 *************************************************************************/

#include <vector>
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "matops.hpp"
#include "vecset.hpp"
#include "v2vtransform.hpp"
#include "transforms.hpp"

using namespace std;

std::vector< ColVector<double> > MotionModel::operator()(const ColVector<double>& x) const {
	std::vector< ColVector<double> > y(1); y[0].resize( x.size() );
	return y;
} // operator()

std::vector< ColVector<double> > SensorModel::operator()(const ColVector<double>& x) const {
	std::vector< ColVector<double> > y(1); y[0].resize( x.size() );
	return y;
} // operator()
