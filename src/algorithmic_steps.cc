/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs the required tracking algorithms.
 *************************************************************************/
#include "algorithmic_steps.hpp"

void predictStates() {
	// Time step and current number of targets
	unsigned N = stateNodes.size();
	unsigned M = currentStates[N-1].size();

	// Resize for indices access
	stateNodes[N].resize(M);
	currentStates[N].resize(M);
	virtualMeasurementVars.resize(M);

	for (unsigned i = 0; i < M; i++) {
		// Get the marginal over previous variables
		currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
		rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);

		// Create a new factor over current variables
		rcptr<Factor> stateJoint = uniqptr<Factor>(new CGM( prevMarginal, 
					initialiseMotionModel(), 
					elementsOfX[currentStates[N][i]],
					initialiseRCovMat()  ));
		stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint));

		// Link Node to preceding node
		stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N-1][i]]);
		stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOfX[currentStates[N-1][i]], prevMarginal);

		// Assign new virtual measurement nodes
		virtualMeasurementVars[i] = addVariables(variables, vecZ, elementsOfZ, mht::kMeasSpaceDim);

		// Create measurement distributions for each measurement space
		predMeasurements[i].resize(mht::kNumSensors); validationRegion[i].resize(mht::kNumSensors);
		for (unsigned j = 0; j < mht::kNumSensors; j++) {
			rcptr<Factor> curMarginal = stateJoint->marginalize(elementsOfX[currentStates[N][i]]);
			curMarginal->inplaceNormalize();

			predMeasurements[i][j] = uniqptr<Factor>(new CGM( curMarginal, 
						mht::kMeasurementModel[j], 
						elementsOfZ[virtualMeasurementVars[i]],
						mht::kQCovMat ) );

			rcptr<Factor> measMarginal = (predMeasurements[i][j])->marginalize(elementsOfZ[virtualMeasurementVars[i]]); 
			validationRegion[i][j]  = (std::dynamic_pointer_cast<CGM>(measMarginal))->momentMatch();
		} // for
	} // for
} // predictStates()

void measurementUpdate() {
	// Time step and current number of targets
	unsigned N = stateNodes.size() - 1;
	unsigned M = currentStates[N].size();
	
	// For every target
	for (unsigned i = 0; i < M; i++) {
		
		// For every sensor
		for (unsigned j = 0; j < mht::kNumSensors; j++) {

			std::vector<ColVector<double>> measurements = measurementManager->getSensorPoints(j, N + 6);
			for (unsigned k = 0; k < measurements.size(); k++) {
				if (measurements[k].size()) {
					double distance = (std::dynamic_pointer_cast<GC>( validationRegion[i][j] ))->mahanalobis(measurements[k]);
					std::cout << distance << std::endl;
				} // if
			} // for

		} // for

	} // for

	virtualMeasurementVars.clear();
	validationRegion.clear();
} // measurementUpdate()
