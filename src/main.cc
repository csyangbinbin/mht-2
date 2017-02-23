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

	initialiseVariables();

	// Step 1 : Get the measurements
	manager = uniqptr<MeasurementManager>(new MeasurementManager("data/test_case_6", mht::kNumSensors));
	kNumberOfTimeSteps = manager->getNumberOfTimeSteps();

	std::vector<ColVector<double>> sample = manager->getSensorPoints(5, 0);

	// Step 2 : Set up the prior
	currentStates[0].clear(); currentStates[0].resize(1);
	currentStates[0][0] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
	
	std::vector<rcptr<Factor>> gc(3);
	gc[0] = uniqptr<Factor>(new GC(elementsOfX[currentStates[0][0]],
				mht::kLaunchStateMean[0],
				mht::kLaunchStateCov[0] ) );


	gc[1] = uniqptr<Factor>(new GC(elementsOfX[currentStates[0][0]],
				mht::kLaunchStateMean[1],
				mht::kLaunchStateCov[1] ) );


	gc[2] = uniqptr<Factor>(new GC(elementsOfX[currentStates[0][0]],
				mht::kLaunchStateMean[2],
				mht::kLaunchStateCov[2] ) );

	rcptr<Factor> prior = uniqptr<Factor>(new CGM(elementsOfX[currentStates[0][0]], gc));
	prior->inplaceNormalize();

	stateNodes[0].clear(); stateNodes[0].resize(1);
	stateNodes[0][0] = uniqptr<Node> (new Node(prior) );

	// Step 3: Loop every time step
	for (unsigned i = 5; i < 6; i++) {
		// Prediction
		predictStates();
		
		// Measurement update
		//measurementUpdate();
		
		// Backwards pass
		
		// Decision making

		// Forwards pass

		// State extraction
	}

	// Step 4: Error metrics
	return 0;
}
