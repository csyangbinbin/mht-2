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
					initialiseMotionModel(), 
					elementsOfX[currentStates[N][i]],
					initialiseRCovMat()  ));
		stateJoint->inplaceNormalize();

		stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint));

		// Link Node to preceding node
		stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N-1][i]]);

		// Assign new virtual measurement nodes
		virtualMeasurementVars[i] = addVariables(variables, vecZ, elementsOfZ, mht::kMeasSpaceDim);

		// Create measurement distributions for each measurement space
		predMeasurements[i].resize(mht::kNumSensors); validationRegion[i].resize(mht::kNumSensors);
		for (unsigned j = 0; j < mht::kNumSensors; j++) {
			std::cout << "SENSOR" << j << std::endl;
			rcptr<Factor> curMarginal = stateJoint->marginalize(elementsOfX[currentStates[N][i]]);
			curMarginal->inplaceNormalize();

			predMeasurements[i][j] = uniqptr<Factor>(new CGM( curMarginal, 
						mht::kMeasurementModel[j], 
						elementsOfZ[virtualMeasurementVars[i]],
						mht::kQCovMat ) );
			
			rcptr<Factor> measMarginal = (predMeasurements[i][j])->marginalize(elementsOfZ[virtualMeasurementVars[i]]); 
			measMarginal->inplaceNormalize();
			
			rcptr<CGM> cgm = (std::dynamic_pointer_cast<CGM>(measMarginal));
			
			validationRegion[i][j]  = cgm->momentMatch();
			std::cout << *validationRegion[i][j] << std::endl;
		}
	}

} // predictStates()

void measurementUpdate() {
} // measurementUpdate()
