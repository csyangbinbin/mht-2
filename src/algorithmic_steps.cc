/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs the required tracking algorithms.
 *************************************************************************/
#include "algorithmic_steps.hpp"

void predictStates() {
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
					mht::kMotionModel, 
					elementsOfX[currentStates[N][i]],
					initialiseRCovMat()  ));
		stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint));

		// Link Node to preceding node
		stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N-1][i]]);
		stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOfX[currentStates[N-1][i]]);

		// Assign new virtual measurement nodes
		rcptr<Factor> curMarginal = stateJoint->marginalize( elementsOfX[currentStates[N][i]] );
		virtualMeasurementVars[i] = addVariables(variables, vecZ, elementsOfZ, mht::kMeasSpaceDim);

		// Create measurement distributions for each measurement space
		predMeasurements[i].resize(mht::kNumSensors); validationRegion[i].resize(mht::kNumSensors);
		for (unsigned j = 0; j < mht::kNumSensors; j++) {
			rcptr<Factor> curMarginal = stateJoint->marginalize(elementsOfX[currentStates[N][i]]);
			predMeasurements[i][j] = uniqptr<Factor>(new CGM( curMarginal, 
						mht::kMeasurementModel[j], 
						elementsOfZ[virtualMeasurementVars[i]],
						mht::kQCovMat ) );

			rcptr<Factor> measMarginal = (predMeasurements[i][j])->marginalize(elementsOfZ[virtualMeasurementVars[i]]); 
			rcptr<CGM> cgm = (std::dynamic_pointer_cast<CGM>(measMarginal));
			validationRegion[i][j]  = cgm->momentMatch();
		}
	}
} // predictStates()

void measurementUpdate() {
	unsigned M = validationRegion[0].size();

	std::cout << M << std::endl;

	for (unsigned i = 0; i < M; i++) {
		std::cout << *validationRegion[0][i] << std::endl;
	}
} // measurementUpdate()
