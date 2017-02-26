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

	// Step 0: Initialise the variables
	initialiseVariables();

	// Step 1 : Get the measurements
	measurementManager = uniqptr<MeasurementManager>(new MeasurementManager("data/test_case_6", mht::kNumSensors));
	kNumberOfTimeSteps = measurementManager->getNumberOfTimeSteps();

	// Step 2 : Create a GraphBuilder object
	graphBuilder = uniqptr<GraphBuilder>(new GraphBuilder());

	// Step 3 : Set up the prior
	currentStates[0].clear(); currentStates[0].resize(1);
	currentStates[0][0] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);

	rcptr<Factor> prior = uniqptr<Factor>(new CGM(elementsOfX[currentStates[0][0]], 
				mht::kGenericWeight,
				mht::kLaunchStateMean,
				mht::kLaunchStateCov));

	rcptr<Node> node = rcptr<Node> (new Node(prior));

	stateNodes[0].clear(); stateNodes[0].resize(1);
	stateNodes[0][0] = uniqptr<Node> (new Node(prior) );

	// Step 4: Loop through every time step
	for (unsigned i = 5; i < 6; i++) {
		// Prediction
		predictStates();
		
		// Measurement update
		measurementUpdate();
		
		// Backwards pass
		
		// Decision making

		// Forwards pass

		// State extraction
	}

	// Step 4: Error metrics
	return 0;
}
