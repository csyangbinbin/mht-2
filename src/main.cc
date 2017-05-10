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
	measurementManager = uniqptr<MeasurementManager>(new MeasurementManager("data/test_case_7", mht::kNumSensors));
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

	// Step 4: Loop through every time step
	for (unsigned i = 1; i < kNumberOfTimeSteps; i++) {
		// Prediction
		predictStates(i, 
				currentStates, 
				virtualMeasurementVars, 
				stateNodes, 
				predMarginals, 
				predMeasurements, 
				validationRegion);

		// Create measurement distributions
		createMeasurementDistributions(i, 
				currentStates, 
				virtualMeasurementVars, 
				stateNodes, 
				measurementNodes, 
				predMarginals, 
				predMeasurements, 
				validationRegion);
		
		// Measurement update
		measurementUpdate(i, 
				stateNodes,
				measurementNodes
				);

		// Backward pass and recalibration
		smoothTrajectory(i, stateNodes);

		// Decision making
		modelSelection(i,
				currentStates, 
				virtualMeasurementVars, 
				stateNodes, 
				measurementNodes, 
				predMarginals, 
				predMeasurements, 
				validationRegion);

		// Forwards pass
		forwardPass(i, stateNodes);

		// Remove states
		removeStates(i, currentStates, stateNodes);
	}

	// State Extraction
	std::cout << "N;x;y;z" << std::endl;
	for (unsigned i = 0; i < kNumberOfTimeSteps; i++) {
		extractStates(i, currentStates, stateNodes);
	} // for

	return 0;
}
