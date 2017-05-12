/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs the required tracking algorithms.
 *************************************************************************/
#include "algorithmic_steps.hpp"

void predictStates(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion
		) {
	// Time step and current number of targets
	unsigned M = stateNodes[N-1].size();

	// Resize for indices access
	stateNodes[N].resize(M); stateNodes[N][0] = 0;
	currentStates[N].resize(M);
	virtualMeasurementVars.resize(M);
	predMarginals.resize(M);

	for (unsigned i = 0; i < M; i++) {
		rcptr<Factor> stateJoint;

		if (i == 0) {
			// Assign new variables to the current state
			currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);

			// Create a clutter state - always has identity of zero
			stateJoint = uniqptr<Factor>(new CGM(elementsOfX[currentStates[N][0]], 
					{1.0},
					{mht::kClutterMean},
					{mht::kClutterCov}));
			stateNodes[N][i] = uniqptr<Node>( new Node(stateJoint, i) );
		} else {
			if (stateNodes[N-1][i] == 0) continue;

			// Assign new variables to the current state
			currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);

			// Get the marginal over previous variables
			rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);
			prevMarginal = std::dynamic_pointer_cast<CGM>(prevMarginal)->momentMatchCGM();

			// Create a new factor over current variables
			stateJoint = uniqptr<Factor>(new CGM( prevMarginal, 
						mht::kMotionModel, 
						elementsOfX[currentStates[N][i]], 
						mht::kRCovMat  ));
			stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint, stateNodes[N-1][i]->getIdentity() ) );

			// Link Node to preceding node
			stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N-1][i]]);
			stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOfX[currentStates[N-1][i]], prevMarginal);
		}

		// Determine the marginal
		predMarginals[i] = stateJoint->marginalize(elementsOfX[currentStates[N][i]]);

		// Assign new virtual measurement nodes
		virtualMeasurementVars[i] = addVariables(variables, vecZ, elementsOfZ, mht::kMeasSpaceDim);

		// Create measurement distributions for each measurement space
		predMeasurements[i].resize(mht::kNumSensors); validationRegion[i].resize(mht::kNumSensors);
		for (unsigned j = 0; j < mht::kNumSensors; j++) {
			predMeasurements[i][j] = uniqptr<Factor>(new CGM( predMarginals[i], 
						mht::kMeasurementModel[j], 
						elementsOfZ[virtualMeasurementVars[i]],
						mht::kQCovMat ) );

			rcptr<Factor> measMarginal = (predMeasurements[i][j])->marginalize(elementsOfZ[virtualMeasurementVars[i]]);
			validationRegion[i][j]  = (std::dynamic_pointer_cast<CGM>(measMarginal))->momentMatch();
		} // for
	} // for
} // predictStates()

void createMeasurementDistributions(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion
		) {
	// Time step and current number of targets
	unsigned M = currentStates[N].size();

	// Clear some things
	measurementNodes[N].clear(); 
	currentMeasurements[N].clear();

	// For each sensor
	for (unsigned i = 0; i < mht::kNumSensors; i++) {
		// Retrieve measurements
		std::vector<ColVector<double>> measurements = measurementManager->getSensorPoints(i, N); // Hack time offset.

		// Local Scope
		std::vector<emdw::RVIdType> sensorMeasurements; 
		std::vector<emdw::RVIdType> associations;

		// Local maps
		std::map<emdw::RVIdType, rcptr<DASS>> assocHypotheses; assocHypotheses.clear();
		std::map<emdw::RVIdType, ColVector<double>> colMeasurements; colMeasurements.clear();

		// For each measurement
		for (unsigned j = 0; j < measurements.size(); j++) { 
			if (measurements[j].size()) {

				// Assign new identity to the measurement
				emdw::RVIdType a = variables.size(); variables.push_back(a);
				emdw::RVIdType z = addVariables(variables, vecZ, elementsOfZ, mht::kMeasSpaceDim);

				// Update local and global scopes
				associations.push_back(a);
				sensorMeasurements.push_back(z);
				currentMeasurements[N].push_back(z);

				assocHypotheses[a] = uniqptr<DASS>(new DASS{0});
				colMeasurements[z] = measurements[j];

				// Form hypotheses over each measurement
				for (unsigned k = 1; k < M; k++) {
					if (validationRegion[k].size() == 0) continue;
					double distance = 
						(std::dynamic_pointer_cast<GC>(validationRegion[k][i]))->mahalanobis(measurements[j]);
						
					if (distance < mht::kValidationThreshold) {
						assocHypotheses[a]->push_back(k); 
					} // if
				} // for
			} // if
		} // for

		if (assocHypotheses.size()) {
			// Determine the 'prior' over association hypotheses
			std::map<emdw::RVIdType, rcptr<Factor>> distributions = graphBuilder->getMarginals(assocHypotheses);
	
			for (unsigned j = 0; j < sensorMeasurements.size(); j++) {
				emdw::RVIdType a = associations[j];
				emdw::RVIdType z = sensorMeasurements[j];

				DASS domain = *assocHypotheses[a];
				unsigned domSize = domain.size();

				//std::cout << "domains: " << domain << std::endl;
				if (domSize > 0) {
					// Create a conditional map for the ConditionalGaussian
					std::map<emdw::RVIdType, rcptr<Factor>> conditionalList; conditionalList.clear();

					for (unsigned k = 0; k < domSize; k++) {
						unsigned p = domain[k];

						// Create a new scope
						emdw::RVIds newScope = predMarginals[p]->getVars();
						newScope.insert(newScope.end(), elementsOfZ[z].begin(), elementsOfZ[z].end());

						// Product of predicted marginals
						conditionalList[p] = uniqptr<Factor>( predMeasurements[p][i]->copy(newScope, false) );

						// Multiply the predicted states by the likelihood function
						for (unsigned l = 0; l < domSize; l++) {
							unsigned q = domain[l];
							if ( q != p ) conditionalList[p]->inplaceAbsorb(predMarginals[q]);
						} // for

						// Introduce evidence into CGM
						conditionalList[p] = conditionalList[p]->observeAndReduce(
								elementsOfZ[z],
								emdw::RVVals{colMeasurements[z][0], colMeasurements[z][1]},
								true);
					} // for
			

					if (N == 13) {
						std::cout << *distributions[a] << std::endl;
					}

					// Create ConditionalGauss - observeAndReduce does work, but for scope reasons this is easier.
					rcptr<Factor> clg = uniqptr<Factor>(new CLG(distributions[a], conditionalList));

					// Create a measurement node and connect it to state nodes
					rcptr<Node> measNode = uniqptr<Node>(new Node(clg));
					for (unsigned k = 0; k < domain.size(); k++) {
						unsigned p = domain[k];
						measNode->addEdge( stateNodes[N][p], predMarginals[p]->getVars(), predMarginals[p]);
						stateNodes[N][p]->addEdge(measNode, predMarginals[p]->getVars());
					} // for
					measurementNodes[N].push_back( measNode );
				} // if
			} // for
		} // if
	} // for

	// Clear temporary values
	predMarginals.clear();
	virtualMeasurementVars.clear();
	validationRegion.clear();
} // createMeasurementDistributions()

void measurementUpdate(const unsigned N,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes
		) {
	//std::cout << "measurementUpdate(" << N << ")" << std::endl;
	unsigned M = measurementNodes[N].size();

	for (unsigned i = 0; i < M; i++) {
		std::vector<std::weak_ptr<Node>> adjacent = measurementNodes[N][i]->getAdjacentNodes();

		for (unsigned j = 0; j < adjacent.size(); j++) {
			// Get the neighbouring state node and message it sent to the measurement clique
			rcptr<Node> stateNode = adjacent[j].lock();
			rcptr<Factor> receivedMessage = measurementNodes[N][i]->getReceivedMessage( stateNode );
			
			// Determine the outgoing message
			emdw::RVIds sepset = measurementNodes[N][i]->getSepset( stateNode );
			rcptr<Factor> outgoingMessage =  (measurementNodes[N][i]->marginalize( sepset, true));

			/*
			if (N == 13) std::cout << stateNode->getIdentity() << std::endl;

			if (stateNode->getIdentity() == 0 && N == 13) {
				double mass = std::dynamic_pointer_cast<CGM>(outgoingMessage)->getMass();
				if (mass > 10) {
					std::cout << *outgoingMessage << std::endl;
				}
			}
			*/

			// Update the factor
			rcptr<Factor> factor = stateNode->getFactor();
			factor->inplaceAbsorb(outgoingMessage);
			factor->inplaceCancel(receivedMessage);

			// Prune and/or merge factor
			std::dynamic_pointer_cast<CGM>(factor)->pruneAndMerge();

			// Update the factor
			stateNode->setFactor(factor);
		} // for
	} // for
} // measurementUpdate()

void smoothTrajectory(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	if (N > mht::kNumberOfBackSteps) { 
		unsigned M = stateNodes[N].size();

		for (unsigned i = 1; i < M; i++) {
			if (stateNodes[N][i] == 0) continue;
			
			// Backwards pass
			for (unsigned j = 0; j < mht::kNumberOfBackSteps; j++) {
				
				// Get the sepset
				emdw::RVIds sepset = stateNodes[N - j][i]->getSepset( stateNodes[N - (j+1)][i] );
			
				// Determine the outgoing message using BUP
				rcptr<Factor> receivedMessage = stateNodes[N-j][i]->getReceivedMessage( stateNodes[N-(j+1)][i] );
				rcptr<Factor> outgoingMessage = stateNodes[N-j][i]->marginalize(sepset, true);
				rcptr<Factor> matched = (std::dynamic_pointer_cast<CGM>(outgoingMessage)->momentMatchCGM());
				
				matched->inplaceCancel(receivedMessage);

				// Received node absorbs message, divides previous message and logs the new message
				stateNodes[N-(j+1)][i]->inplaceAbsorb( matched.get()  );
				stateNodes[N-(j+1)][i]->logMessage( stateNodes[N-j][i], uniqptr<Factor>(matched->copy()) );
			} // for
		} // for
	} // if
} // smoothTrajectory()

void modelSelection(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion
		) {	
	if ( N > mht::kNumberOfBackSteps + 1 ) {
		//std::cout << "modelSelection()" << std::endl;

		unsigned K = N - mht::kNumberOfBackSteps;
		unsigned M = stateNodes[K-1].size();

		// Number of active targets
		unsigned numberOfTargets = 0;
		for (unsigned i = 1; i < M; i++) {
			if (stateNodes[K][i] != 0) numberOfTargets += 1;
		} // for

		// Determine odds for current model
		if (N == 15)std::cout << "ModelOneOdds" << std::endl;
		double modelOneOdds = calculateEvidence(K, stateNodes);		

		// Create new model - just copies of stateNodes and newMeasurementNodes
		std::map<unsigned, emdw::RVIds> newCurrentStates; newCurrentStates[K-1].resize(M);
		std::map<unsigned, std::vector<rcptr<Node>>> newStateNodes; newStateNodes[K-1].resize(M);
		std::map<unsigned, std::vector<rcptr<Node>>> newMeasurementNodes;

		//std::cout << "New prior creation: " << std::endl;
		
		// Copy the stateNodes for time K-1
		for (unsigned i = 0; i < M; i++) {
			newCurrentStates[K-1][i] = currentStates[K-1][i];
			newStateNodes[K-1][i] = stateNodes[K-1][i];
		} // for
			
		// Create a prior for the new target and add it to the preceding time step
		newCurrentStates[K-1].push_back(addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim));
		rcptr<Factor> newTargetPrior = uniqptr<Factor>(new CGM(elementsOfX[newCurrentStates[K-1][M]], 
					{mht::kGenericWeight[0]},
					{mht::kLaunchStateMean[0]},
					{mht::kLaunchStateCov[0]}));
		newStateNodes[K-1].push_back(uniqptr<Node> (new Node(newTargetPrior, M) ));

		// Propagate the new model forward
		for (unsigned i = K; i <= N; i++) {
			// Predict states
			predictStates(i, newCurrentStates, virtualMeasurementVars, newStateNodes, predMarginals, predMeasurements, 
				validationRegion);

			// Recreate measurement distributions
			createMeasurementDistributions(i, newCurrentStates, virtualMeasurementVars, newStateNodes,
					newMeasurementNodes, predMarginals, predMeasurements, validationRegion);

			// Measurement update
			measurementUpdate(i, newStateNodes, newMeasurementNodes);
	
			/*
			if (i == 13 && N == 15) {
				std::cout << "i: " << i << std::endl;
				std::cout << *newStateNodes[i][0] << std::endl;
			}
			*/

		} // for
		//smoothTrajectory(N, newStateNodes);

		// Calculate new model odds
		if (N==15) std::cout << "ModelTwoOdds" << std::endl;
		double modelTwoOdds = calculateEvidence(K, newStateNodes) + log(mht::kTimeStep*3) - log(numberOfTargets+1);
		//std::cout << "modelTwoOdds: " << exp(modelTwoOdds) << std::endl;
		
		if (modelTwoOdds > modelOneOdds) {
			std::cout << "K: " << K << "- Use model two now!" << std::endl;	
			std::cout << "modelOneOdds: " << modelOneOdds << std::endl;
			std::cout << "modelTwoOdds: " << modelTwoOdds << std::endl;

			/*
			// Replace model one
			for (unsigned i = K-1; i <= N; i++) {

				unsigned M = newStateNodes[i].size();

				stateNodes[i].clear(); stateNodes[i].resize(M);
				currentStates[i].clear(); currentStates[i].resize(M);

				for (unsigned j = 0; j < M; j++) {
					stateNodes[i][j] = newStateNodes[i][j];
					currentStates[i][j] = newCurrentStates[i][j];
				} // for
			} // for
			*/
		} // if

		// Clear model two
		// newStateNodes.clear(); newCurrentStates.clear();
		
		//std::cout << "\n\n" << std::endl;
	} // if
} // modelSelection()

void forwardPass(unsigned const N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	if (N > mht::kNumberOfBackSteps) {
		unsigned M = stateNodes[N].size();

		for (unsigned i = 1; i < M; i++) {
			for (unsigned j = mht::kNumberOfBackSteps; j > 0; j--) {
				if (stateNodes[N][i] == 0) continue;

				// Get current sepsets and scope
				emdw::RVIds presentVars = stateNodes[N-j][i]->getSepset( stateNodes[N - (j-1)][i] );

				// Determine the outgoing message using BUP
				rcptr<Factor> receivedMessage = stateNodes[N-j][i]->getReceivedMessage( stateNodes[N-(j-1)][i] );
				rcptr<Factor> outgoingMessage = stateNodes[N-j][i]->marginalize(presentVars, true);
				
				rcptr<Factor> matched = (std::dynamic_pointer_cast<CGM>(outgoingMessage)->momentMatchCGM());
				matched->inplaceCancel(receivedMessage);

				// Receiving node absorbs and logs the message
				stateNodes[N-(j-1)][i]->absorb(matched.get() );
				stateNodes[N-(j-1)][i]->logMessage( stateNodes[N-j][i], matched );
			} // for	
		} // for
	} // if
} // forwardPass()


void removeStates(const unsigned N, 
		std::map<unsigned, emdw::RVIds>& currentStates,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes
		) {
	unsigned M = stateNodes[N].size();

	for (unsigned i = 1; i < M; i++) {
		if (stateNodes[N][i] == 0) continue;
		rcptr<Factor> marginal = stateNodes[N][i]->marginalize(elementsOfX[currentStates[N][i]], true);
		rcptr<Factor> matched = std::dynamic_pointer_cast<CGM>( marginal )->momentMatch();
	
		ColVector<double> mean =  std::dynamic_pointer_cast<GC>(matched)->getMean();
		if (mean[4] < -7.0) stateNodes[N][i] = 0;
	} // for

} // removeStates()
