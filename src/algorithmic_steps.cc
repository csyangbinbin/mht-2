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
	
	stateNodes[N].resize(M);
	currentStates[N].resize(M);
	virtualMeasurements.resize(kNumSensors);

	for (unsigned i = 0; i < M; i++) {
		// Get the marginal over previous variables
		currentStates[N][i] = addVariables(variables, vecX, elementsOfX, kStateSpaceDim);
		rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);

		// Create a new factor over current variables
		rcptr<Factor> stateJoint = uniqptr<Factor>(new CGM( prevMarginal, kMotionModel, elementsOfX[currentStates[N][i]], kRCovMat  ));
		stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint));

		// Link Node to preceding node
		stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N][i]]);
		stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOfX[currentStates[N][i]]);

		// Create measurement distributions for each measurement space
		rcptr<Factor> curMarginal = stateJoint->marginalize( elementsOfX[currentStates[N][i]] );
		virtualMeasurements[i].push_back(addVariables(variables, vecZ, elementsOfZ, kMeasSpaceDim));
		
		predMeasurements[i].resize(kNumSensors); validationRegion[i].resize(kNumSensors);
		for (unsigned j = 0; j < kNumSensors; j++) {
			predMeasurements[i][j] = uniqptr<Factor>(new CGM( curMarginal, kMeasurementModel[j], elementsOfZ[virtualMeasurements[i][0]], kQCovMat ));
			rcptr<Factor> measMarginal = predMeasurements[i][j]->marginalize(elementsOfZ[virtualMeasurements[i][0]]); 
			validationRegion[i][j] = (std::dynamic_pointer_cast<CGM>( measMarginal ))->momentMatch();
		}
	}

	virtualMeasurements.clear();

} // predictStates()
