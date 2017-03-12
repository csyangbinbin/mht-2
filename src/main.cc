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
	currentStates[0].clear(); currentStates[0].resize(2); vecX.push_back(0);
	stateNodes[0].clear(); stateNodes[0].resize(2); stateNodes[0][0] = 0;

	// Tee 1
	currentStates[0][1] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
	rcptr<Factor> teeOne = uniqptr<Factor>(new CGM(elementsOfX[currentStates[0][1]], 
				{mht::kGenericWeight[0]},
				{mht::kLaunchStateMean[0]},
				{mht::kLaunchStateCov[0]}));
	stateNodes[0][1] = uniqptr<Node> (new Node(teeOne, 1) );

	// Tee 2
	/*
	currentStates[0][2] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
	rcptr<Factor> teeTwo = uniqptr<Factor>(new CGM(elementsOfX[currentStates[0][2]], 
				{mht::kGenericWeight[1]},
				{mht::kLaunchStateMean[1]},
				{mht::kLaunchStateCov[1]}));
	stateNodes[0][2] = uniqptr<Node> (new Node(teeTwo, 2) );

	// Tee 3
	currentStates[0][3] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
	rcptr<Factor> teeThree = uniqptr<Factor>(new CGM(elementsOfX[currentStates[0][3]], 
				{mht::kGenericWeight[2]},
				{mht::kLaunchStateMean[2]},
				{mht::kLaunchStateCov[2]}));
	stateNodes[0][3] = uniqptr<Node> (new Node(teeThree, 3) );
	*/

	std::cout << "N;x;y;z" << std::endl;

	// Step 4: Loop through every time step
	for (unsigned i = 1; i <= 145; i++) {
		// Prediction
		predictStates(i);
		
		// Create measurement distributions
		createMeasurementDistributions(i);
		
		// Measurement update
		measurementUpdate(i);
		
		// Backward pass and recalibration
		smoothTrajectory(i);

		// Decision making
		
		// Forwards pass

		// State extraction
	}

	// Step 4: Error metrics
	return 0;
}
