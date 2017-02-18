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

	for (unsigned i = 0; i < M; i++) {
		currentStates[N][i] = addVariables(variables, vecX, elementsOf, kStateSpaceDim);
		rcptr<Factor> marginal = (stateNodes[N-1][i])->marginalize(elementsOf[currentStates[N-1][i]]);

		rcptr<Factor> joint = uniqptr<Factor>(new CGM( marginal, kMotionModel, elementsOf[currentStates[N][i]], kRCovMat  ));
		stateNodes[N][i] = uniqptr<Node> (new Node(joint));
		
		stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOf[currentStates[N][i]]);
		stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOf[currentStates[N][i]]);
	}

} // predictStates()
