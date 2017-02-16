/*************************************************************************
 *  Compilation: ./run_main 0
 *  Execution: ./run_main 0
 *  Dependencies:
 *
 * Main app, runs everything.
 *************************************************************************/

#include <iostream>
#include <map>
#include "genvec.hpp"
#include "genmat.hpp"
#include "anytype.hpp"
#include "emdw.hpp"
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "linear_gaussian.hpp"
#include "transforms.hpp"
#include "utils.hpp"
#include "system_constants.hpp"

/**
 * Main app, runs small examples for now.
 *
 * @author SCJ Robertson
 * @since 03/10/16
 */
int main(int, char *argv[]) {
	
	emdw::RVIds currentX;
	currentX.push_back(addVariables(variables, vecX, elementsOf, N));

	rcptr<Factor> prior = uniqptr<Factor>(new CanonicalGaussianMixture(elementsOf[currentX[0]], {1.0}, {kLaunchStateMean[0]}, {kLaunchStateCov[0]}));

	return 0;
}
