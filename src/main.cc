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
	unsigned operationMode = 1;

	// Step 1 : Get the measurements
	measurementManager = uniqptr<MeasurementManager>(new MeasurementManager("data/test_case_6", mht::kNumSensors));
	kNumberOfTimeSteps = measurementManager->getNumberOfTimeSteps();

	// Step 2 : Create a GraphBuilder object
	graphBuilder = uniqptr<GraphBuilder>(new GraphBuilder());

	// Step 3 : Set up the prior
	currentStates[0].clear(); currentStates[0].resize(mht::kNumSensors+1); 
	stateNodes[0].clear(); stateNodes[0].resize(mht::kNumSensors+1); 
	
	for (unsigned i = 0; i < mht::kNumSensors; i++) { 
		vecX.push_back(i);
		stateNodes[0][i] = 0;
	} // for

	// Tee 1
	currentStates[0][mht::kNumSensors] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
	rcptr<Factor> teeOne = uniqptr<Factor>(new CGM(elementsOfX[currentStates[0][mht::kNumSensors]], 
				{1.0},
				{1.0*mht::kGenericMean},
				{1.0*mht::kGenericCov}));
	stateNodes[0][mht::kNumSensors] = uniqptr<Node> (new Node(teeOne, mht::kNumSensors) );

	// Step 4: Loop through every time step
	for (unsigned i = 1; i < kNumberOfTimeSteps; i++) {

		if (operationMode == 0) {
			// Prediction
			predictStatesSU(i, 
					currentStates, 
					virtualMeasurementVars, 
					stateNodes, 
					predMarginals, 
					predMeasurements, 
					validationRegion);

			// Create measurement distributions
			createMeasurementDistributionsSU(i, 
					currentStates, 
					virtualMeasurementVars, 
					stateNodes, 
					measurementNodes, 
					predMarginals, 
					predMeasurements, 
					validationRegion);
			
			// Measurement update
			measurementUpdateSU(i, 
					stateNodes,
					measurementNodes
					);

			// Backward pass and recalibration
			smoothTrajectory(i, stateNodes);

			// Decision making
			modelSelectionSU(i,
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

		} else if (operationMode == 1) {
			// Prediction
			predictStatesAU(i,
					currentStates, 
					stateNodes);

			// Measurement update
			measurementUpdateAU(i, 
					currentStates, 
					virtualMeasurementVars, 
					stateNodes, 
					measurementNodes, 
					predMarginals, 
					predMeasurements, 
					validationRegion);

			// Backward pass and recalibration
			smoothTrajectory(i, stateNodes);

			// Decision making
			modelSelectionAU(i,
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
		} // if
	}

	// State Extraction
	for (unsigned i = 0; i < kNumberOfTimeSteps; i++) {
		extractStates(i, currentStates, stateNodes);
	} // for

	return 0;
}
