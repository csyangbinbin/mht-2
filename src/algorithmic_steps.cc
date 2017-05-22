/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs the required tracking algorithms.
 *************************************************************************/
#include "algorithmic_steps.hpp"

void predictStatesSU(const unsigned N,
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
	stateNodes[N].resize(M); 
	currentStates[N].resize(M);
	virtualMeasurementVars.resize(M);
	predMarginals.resize(M);

	for (unsigned i = 0; i < M; i++) {
		rcptr<Factor> stateJoint;

		if (i < mht::kNumSensors) {
			// Assign new variables to the current state
			currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);

			// Create a clutter state - always has identity of a sensor
			stateJoint = uniqptr<Factor>(new CGM(elementsOfX[currentStates[N][i]], 
					{1.0},
					{mht::kClutterMean[i]},
					{mht::kClutterCov[i]}));
			stateNodes[N][i] = uniqptr<Node>( new Node(stateJoint, i) );
		} else {
			if (stateNodes[N-1][i] == nullptr) continue;

			// Assign new variables to the current state
			currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);

			// Get the marginal over previous variables
			rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);
			prevMarginal->inplaceNormalize();
			
			// Moment match the distribution for the received message
			rcptr<Factor> receivedMessage = std::dynamic_pointer_cast<CGM>(prevMarginal)->momentMatchCGM();
			
			// Create a new factor over current variables
			stateJoint = uniqptr<Factor>(new CGM( prevMarginal, 
						mht::kMotionModel, 
						elementsOfX[currentStates[N][i]], 
						mht::kRCovMat  ));
			stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint, stateNodes[N-1][i]->getIdentity() ) );

			// Link Node to preceding node
			stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N-1][i]]);
			stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOfX[currentStates[N-1][i]], receivedMessage);
		} // if

		// Determine the marginal
		predMarginals[i] = stateJoint->marginalize(elementsOfX[currentStates[N][i]]);
		predMarginals[i] = std::dynamic_pointer_cast<CGM>(predMarginals[i])->momentMatchCGM();

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
} // predictStatesSU()


void predictStatesAU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes
		) {
	// Time step and current number of targets
	unsigned M = stateNodes[N-1].size();

	// Resize for indices access
	stateNodes[N].resize(M); stateNodes[N][0] = 0;
	currentStates[N].resize(M);

	for (unsigned i = 0; i < M; i++) {
		rcptr<Factor> stateJoint;

		if (i < mht::kNumSensors) {
			// Assign new variables to the current state
			currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);

			// Create a clutter state - always has identity of a sensor
			stateJoint = uniqptr<Factor>(new CGM(elementsOfX[currentStates[N][i]], 
					{1.0},
					{mht::kClutterMean[i]},
					{mht::kClutterCov[i]}));
			stateNodes[N][i] = uniqptr<Node>( new Node(stateJoint, i) );
		} else {
			if (stateNodes[N-1][i] == nullptr) continue;

			// Assign new variables to the current state
			currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);

			// Get the marginal over previous variables
			rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);
			prevMarginal->inplaceNormalize();
			
			// Moment match the distribution for the received message
			rcptr<Factor> receivedMessage = std::dynamic_pointer_cast<CGM>(prevMarginal)->momentMatchCGM();

			// Create a new factor over current variables
			stateJoint = uniqptr<Factor>(new CGM( prevMarginal, 
						mht::kMotionModel, 
						elementsOfX[currentStates[N][i]], 
						mht::kRCovMat  ));
			stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint, stateNodes[N-1][i]->getIdentity() ) );

			// Link Node to preceding node
			stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N-1][i]]);
			stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOfX[currentStates[N-1][i]], receivedMessage);
		} // if
	} // for
} // predictStatesAU()


void createMeasurementDistributionsSU(const unsigned N,
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

	// For each sensorI
	for (unsigned i = 0; i < mht::kNumSensors; i++) {
		// Retrieve measurements
		std::vector<ColVector<double>> measurements = measurementManager->getSensorPoints(i, N + 6); // Hack time offset.

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

				assocHypotheses[a] = uniqptr<DASS>(new DASS{ (unsigned short) i  });
				colMeasurements[z] = measurements[j];

				// Form hypotheses over each measurement
				for (unsigned k = mht::kNumSensors; k < M; k++) {
					if (validationRegion[k].size() == 0) continue;
					
					//std::cout << "k : " << k << std::endl;
					double distance = 
						(std::dynamic_pointer_cast<GC>(validationRegion[k][i]))->mahalanobis(measurements[j]);
					//std::cout << "Distance : " << distance  << std::endl;

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
} // createMeasurementDistributionsSU()


void createMeasurementDistributionsAU(const unsigned N,
		const unsigned sensorNumber,
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
	measurementNodes[N].clear(); currentMeasurements[N].clear();
	virtualMeasurementVars.resize(M); predMarginals.resize(M);

	// Create predicted measurement distributions for a particular sensor
	for (unsigned i = 0; i < M; i++) {
		if (stateNodes[N][i] == nullptr) continue;
		if (i != sensorNumber && i < mht::kNumSensors) continue;

		// Determine the predicted marginal
		predMarginals[i] = stateNodes[N][i]->marginalize(elementsOfX[currentStates[N][i]]);
		predMarginals[i] = std::dynamic_pointer_cast<CGM>(predMarginals[i])->momentMatchCGM();
		
		// Add in place holder values for measurement vars
		virtualMeasurementVars[i] = addVariables(variables, vecZ, elementsOfZ, mht::kMeasSpaceDim);

		// Create the predicted measurement distribution
		predMeasurements[i].resize(1); validationRegion[i].resize(1);
		predMeasurements[i][0] = uniqptr<Factor>(new CGM( predMarginals[i], 
					mht::kMeasurementModel[sensorNumber], 
					elementsOfZ[virtualMeasurementVars[i]],
					mht::kQCovMat ) );

		// Cast the predicted distribution to a GaussCannical so Mahalanobis distance can be used.
		rcptr<Factor> measMarginal = (predMeasurements[i][0])->marginalize(elementsOfZ[virtualMeasurementVars[i]]);
		validationRegion[i][0]  = (std::dynamic_pointer_cast<CGM>(measMarginal))->momentMatch();
	} // for

	// Get the measurement for this particular sensor
	std::vector<ColVector<double>> measurements = measurementManager->getSensorPoints(sensorNumber, N + mht::timeOffSet); // Hack time offset.

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

			assocHypotheses[a] = uniqptr<DASS>(new DASS{ (unsigned short) sensorNumber  });
			colMeasurements[z] = measurements[j];

			// Form hypotheses over each measurement
			for (unsigned k = mht::kNumSensors; k < M; k++) {
				if (validationRegion[k].size() == 0) continue;

				double distance = 
					(std::dynamic_pointer_cast<GC>(validationRegion[k][0]))->mahalanobis(measurements[j]);

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

			if (domSize > 0) {
				// Create a conditional map for the ConditionalGaussian
				std::map<emdw::RVIdType, rcptr<Factor>> conditionalList; conditionalList.clear();

				for (unsigned k = 0; k < domSize; k++) {
					unsigned p = domain[k];

					// Create a new scope
					emdw::RVIds newScope = predMarginals[p]->getVars();
					newScope.insert(newScope.end(), elementsOfZ[z].begin(), elementsOfZ[z].end());

					// Product of predicted marginals
					conditionalList[p] = uniqptr<Factor>( predMeasurements[p][0]->copy(newScope, false) );

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

	// Clear temporary values
	predMarginals.clear();
	virtualMeasurementVars.clear();
	validationRegion.clear();
} // createMeasurementDistributionsAU()

void measurementUpdateSU(const unsigned N,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes
		) {
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
} // measurementUpdateSU()

void measurementUpdateAU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion
		) {
		
	for (unsigned i = 0; i < mht::kNumSensors; i++) {
		createMeasurementDistributionsAU(N,
				i,
				currentStates, 
				virtualMeasurementVars, 
				stateNodes,
				measurementNodes, 
				predMarginals, 
				predMeasurements, 
				validationRegion);

		unsigned M = measurementNodes[N].size();

		for (unsigned j = 0; j < M; j++) {
			std::vector<std::weak_ptr<Node>> adjacent = measurementNodes[N][j]->getAdjacentNodes();

			for (unsigned k = 0; k < adjacent.size(); k++) {
				// Get the neighbouring state node and message it sent to the measurement clique
				rcptr<Node> stateNode = adjacent[k].lock();
				rcptr<Factor> receivedMessage = measurementNodes[N][j]->getReceivedMessage( stateNode );

				// Determine the outgoing message
				emdw::RVIds sepset = measurementNodes[N][j]->getSepset( stateNode );
				rcptr<Factor> outgoingMessage =  (measurementNodes[N][j]->marginalize( sepset, true));

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
	} // for
} // createMeasurementDistributionsAU()

void smoothTrajectory(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	if (N > mht::kNumberOfBackSteps) { 
		unsigned M = stateNodes[N].size();

		for (unsigned i = mht::kNumSensors; i < M; i++) {
			if (stateNodes[N][i] == nullptr) continue;
			
			// Backwards pass
			for (unsigned j = 0; j < mht::kNumberOfBackSteps; j++) {
				
				// Get the sepset
				emdw::RVIds sepset = stateNodes[N - j][i]->getSepset( stateNodes[N - (j+1)][i] );
			
				// Determine the outgoing message using BUP
				rcptr<Factor> receivedMessage = stateNodes[N-j][i]->getReceivedMessage( stateNodes[N-(j+1)][i] );
				rcptr<Factor> outgoingMessage = stateNodes[N-j][i]->marginalize(sepset, true);
				rcptr<Factor> matched = (std::dynamic_pointer_cast<CGM>(outgoingMessage)->momentMatchCGM());
				
				matched->inplaceCancel(receivedMessage);

				stateNodes[N-(j+1)][i]->inplaceAbsorb( matched.get()  );
				stateNodes[N-(j+1)][i]->logMessage( stateNodes[N-j][i], uniqptr<Factor>(matched->copy()) );	
			} // for
		} // for
	} // if
} // smoothTrajectory()

void modelSelectionSU(const unsigned N,
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
		for (unsigned i = mht::kNumSensors; i < M; i++) {
			if (stateNodes[K][i] != nullptr) numberOfTargets += 1;
		} // for

		if (numberOfTargets < mht::maxNumberOfTargets) {   

			// Determine odds for current model
			double modelOneOdds = calculateEvidence(K, K, stateNodes);		

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
						{1.0},
						{1.0*mht::kGenericMean},
						{1.0*mht::kGenericCov}));
			newStateNodes[K-1].push_back(uniqptr<Node> (new Node(newTargetPrior, M) ));

			// Propagate the new model forward
			for (unsigned i = K; i <= N; i++) {
				// Predict states
				predictStatesSU(i, 
						newCurrentStates, 
						virtualMeasurementVars, 
						newStateNodes, 
						predMarginals, 
						predMeasurements, 
						validationRegion);

				// Recreate measurement distributions
				createMeasurementDistributionsSU(i, 
						newCurrentStates, 
						virtualMeasurementVars, 
						newStateNodes,
						newMeasurementNodes, 
						predMarginals, 
						predMeasurements, 
						validationRegion);

				// Measurement update
				measurementUpdateSU(i, 
						newStateNodes, 
						newMeasurementNodes);
			} // for
			smoothTrajectory(N, newStateNodes);

			// Calculate new model odds
			double modelTwoOdds = calculateEvidence(K, K, newStateNodes) 
				+ log(mht::kTimeStep*1) - log(numberOfTargets+1);
			
			if (modelTwoOdds > modelOneOdds) {
				std::cerr << "K: " << K << "- Use model two now!" << "\n";	
				std::cerr << "modelOneOdds: " << modelOneOdds << "\n";
				std::cerr << "modelTwoOdds: " << modelTwoOdds << "\n";

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
			} // if
	   } // if
	} // if
} // modelSelectionSU()

void modelSelectionAU(const unsigned N,
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
		for (unsigned i = mht::kNumSensors; i < M; i++) {
			if (stateNodes[K][i] != nullptr) numberOfTargets += 1;
		} // for

		if (numberOfTargets < mht::maxNumberOfTargets) {    //(K-1) > 9 && (K-1) < 20) {

			// Determine odds for current model
			double modelOneOdds = calculateEvidence(K, N, stateNodes);		

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
						{1.0},
						{1.0*mht::kGenericMean},
						{1.0*mht::kGenericCov}));
			newStateNodes[K-1].push_back(uniqptr<Node> (new Node(newTargetPrior, M) ));

			// Propagate the new model forward
			for (unsigned i = K; i <= N; i++) {
				// Predict states
				predictStatesAU(i, 
						newCurrentStates, 
						newStateNodes);

				// Recreate measurement distributions
				measurementUpdateAU(i, 
						newCurrentStates, 
						virtualMeasurementVars, 
						newStateNodes,
						newMeasurementNodes, 
						predMarginals, 
						predMeasurements, 
						validationRegion);
			} // for
			smoothTrajectory(N, newStateNodes);

			// Calculate new model odds
			double modelTwoOdds = calculateEvidence(K, N, newStateNodes) 
				+ log(mht::kTimeStep*1) - log(numberOfTargets+1);
			
			if (modelTwoOdds > modelOneOdds) {
				std::cerr << "K: " << K << "- Use model two now!" << "\n";	
				std::cerr << "modelOneOdds: " << modelOneOdds << "\n";
				std::cerr << "modelTwoOdds: " << modelTwoOdds << "\n";

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
			} // if
	   } // if
	} // if
} // modelSelectionAU()

void forwardPass(unsigned const N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	if (N > mht::kNumberOfBackSteps) {
		unsigned M = stateNodes[N].size();

		for (unsigned i = mht::kNumSensors; i < M; i++) {
			if (stateNodes[N][i] == nullptr) continue;

			for (unsigned j = mht::kNumberOfBackSteps; j > 0; j--) {
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

	for (unsigned i = mht::kNumSensors; i < M; i++) {
		if (stateNodes[N][i] == nullptr) continue;
		rcptr<Factor> marginal = stateNodes[N][i]->marginalize(elementsOfX[currentStates[N][i]], true);
		rcptr<Factor> matched = std::dynamic_pointer_cast<CGM>( marginal )->momentMatch();
	
		ColVector<double> mean =  std::dynamic_pointer_cast<GC>(matched)->getMean();
		if (mean[4] < -7.0) {
			stateNodes[N][i] = nullptr;
			std::cerr << "N: " << N << ", removed Target " << i << "\n";
		} // if
	} // for

} // removeStates()
