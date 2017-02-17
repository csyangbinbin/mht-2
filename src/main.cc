/*************************************************************************
 *  Compilation: ./run_main 0
 *  Execution: ./run_main 0
 *  Dependencies:
 *
 * Main app, runs everything.
 *************************************************************************/
#include "system_constants.hpp"
#include "algorithmic_steps.hpp"
#include "utils.hpp"

/**
 * Main app, runs small examples for now.
 *
 * @author SCJ Robertson
 * @since 03/10/16
 */
int main(int, char *argv[]) {
	// Step 1 : Get the measurements
	manager = uniqptr<MeasurementManager>(new MeasurementManager("data/test_case_6", kNumSensors));
	kNumberOfTimeSteps = manager->getNumberOfTimeSteps();

	// Step 2 : Set up the prior
	emdw::RVIds currentX;
	currentX.push_back(addVariables(variables, vecX, elementsOf, N));
	rcptr<Factor> prior = uniqptr<Factor>(new CGM(elementsOf[currentX[0]], {1.0}, {kLaunchStateMean[0]}, {kLaunchStateCov[0]}));

	nodes[0].push_back( uniqptr<Node> (new Node(prior) ) );
	factors[0].push_back(prior);

	// Step 3: Loop every time step
	for (unsigned i = 5; i < 6; i++) {

		// Prediction
		predictStates();
		
		// Measurement update
		
		// Backwards pass
		
		// Decision making

		// Forwards pass

		// State extraction
	}

	// Step 4: Error metrics
	return 0;
}
