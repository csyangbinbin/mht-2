/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 I* Runs the required tracking algorithms.
 *************************************************************************/
#include "algorithmic_steps.hpp"

void predictStates(const unsigned N) {
	// Time step and current number of targets
	unsigned M = currentStates[N-1].size();

	// Resize for indices access
	stateNodes[N].resize(M); stateNodes[N][0] = 0;
	currentStates[N].resize(M);
	virtualMeasurementVars.resize(M);
	predMarginals.resize(M);

	for (unsigned i = 0; i < M; i++) {
		// Assign new variables to the current state
		currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
		
		rcptr<Factor> stateJoint;
		if (i == 0) {
			// Create a clutter state - always has identity of zero
			stateJoint = uniqptr<Factor>(new CGM(elementsOfX[currentStates[N][0]], 
					{1.0},
					{mht::kClutterMean},
					{mht::kClutterCov}));
			stateNodes[N][i] = uniqptr<Node>( new Node(stateJoint, i) );
		} else {
			// Get the marginal over previous variables
			rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);
			
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
		//predMarginals[i]->inplaceNormalize();
		
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

void createMeasurementDistributions(const unsigned N) {
	//std::cout << "createMeasurementDistributions()" << std::endl;
	
	// Time step and current number of targets
	unsigned M = currentStates[N].size();
	//std::cout << "M = " << M << std::endl;

	// Clear some things
	measurementNodes[N].clear(); 
	currentMeasurements[N].clear();

	// For each sensor
	for (unsigned i = 0; i < mht::kNumSensors; i++) {
		// Retrieve measurements
		std::vector<ColVector<double>> measurements = measurementManager->getSensorPoints(i, N + 4); // Hack time offset.

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
							if ( q != p ) conditionalList[p] = 
								conditionalList[p]->absorb(predMarginals[q] );
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
} // createMeasurementDistributions()

void measurementUpdate(const unsigned N) {
	unsigned M = measurementNodes[N].size();

	for (unsigned i = 0; i < M; i++) {
		std::vector<std::weak_ptr<Node>> adjacent = measurementNodes[N][i]->getAdjacentNodes();

		for (unsigned j = 0; j < adjacent.size(); j++) {
			// Get the neighbouring state node and message it sent to the measurement clique
			rcptr<Node> stateNode = adjacent[j].lock();
			if (stateNode->getIdentity() != 0) {
				rcptr<Factor> receivedMessage = measurementNodes[N][i]->getReceivedMessage( stateNode );

				// Determine the outgoing message
				emdw::RVIds sepset = measurementNodes[N][i]->getSepset( stateNode );
				rcptr<Factor> outgoingMessage =  measurementNodes[N][i]->marginalize( sepset, true );

				// Get the factor, absorb  the incoming message and then prune
				rcptr<Factor> factor = stateNode->getFactor();
				factor->inplaceAbsorb( outgoingMessage );
				factor->inplaceCancel( receivedMessage );

				//rcptr<Factor> marg = factor->marginalize( sepset, true );
				//std::cout << *factor << std::endl;

				// if (stateNode->getIdentity() != 0) std::dynamic_pointer_cast<CGM>(factor)->pruneAndMerge();
				
				// Update the factor
				stateNode->setFactor(factor);
				stateNode->logMessage( measurementNodes[N][i] ,outgoingMessage );
			}
		} // for
	} // for
} // measurementUpdate()

void smoothTrajectory(const unsigned N) {
	if (N > mht::kNumberOfBackSteps) { 
		unsigned M = stateNodes[N].size();

		for (unsigned i = 1; i < M; i++) {
			// Backwards pass
			for (unsigned j = 0; j < mht::kNumberOfBackSteps; j++) {
				// Get the sepset
				emdw::RVIds sepset = stateNodes[N - j][i]->getSepset( stateNodes[N - (j+1)][i] );
			
				// Determine the outgoing message using BUP
				rcptr<Factor> receivedMessage = stateNodes[N-j][i]->getReceivedMessage( stateNodes[N-(j+1)][i] );
				rcptr<Factor> outgoingMessage = stateNodes[N-j][i]->marginalize(sepset, true);

				// Received node absorbs messgae, divides previous message and logs the new message
				stateNodes[N-(j+1)][i]->inplaceAbsorb( outgoingMessage.get()  );
				stateNodes[N-(j+1)][i]->inplaceCancel( receivedMessage.get() );
				stateNodes[N-(j+1)][i]->logMessage( stateNodes[N-j][i], outgoingMessage  );
			} // for
		} // for
	} // if
} // smoothTrajectory()

void modelSelection(const unsigned N) {	
	if ( N > mht::kNumberOfBackSteps ) {
		// Number of targets a few steps back.
		unsigned K = N - mht::kNumberOfBackSteps;
		unsigned M = stateNodes[K].size();
		
		// Determine the evidence based on the current model
		double modelOneOdds = calculateEvidence(K);
		std::cout << "Current Model's evidence: " << exp(modelOneOdds) << std::endl;

		// Cache the current model
		std::vector<rcptr<Node>> stateNodesCache(M); 
		emdw::RVIds currentStatesCache(M);

		for (unsigned i = 0; i < M; i++) {
			stateNodesCache[i] = stateNodes[K][i];
			currentStatesCache[i] = currentStates[K][i];
		} // for

		// Add in a prior for a new target
		M = stateNodes[K-1].size();
		currentStates[K-1].push_back(addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim));
		rcptr<Factor> prior = uniqptr<Factor>(new CGM(elementsOfX[currentStates[K-1][M]], 
					{1.0},
					{mht::kGenericMean},
					{mht::kGenericCov}));
		stateNodes[K-1].push_back( uniqptr<Node>( new Node(prior, M) ) );

		// Re-predict the states
		predictStates(K);
		createMeasurementDistributions(K);
		measurementUpdate(K);

		// Calculate the current odds
		double modelTwoOdds = calculateEvidence(K);
		std::cout << "New Model's evidence: " << exp(modelTwoOdds) << std::endl;

		// Determine the odds ratio - only choose model if it is strictly more likely.
		if (true) {
			// Clear the state nodes and current state variables
			stateNodes[K].clear(); currentStates[K].clear();
			stateNodes[K].resize(M); currentStates[K].resize(M);
			
			// Reassign the model
			for (unsigned i = 0; i < M; i++) {
				stateNodes[K][i] = stateNodesCache[i];
				currentStates[K][i] = currentStatesCache[i];
			} // for
		} // if 
	} // if
} // modelSelection()

void forwardPass(unsigned const N) {
	if (N > mht::kNumberOfBackSteps) {
		unsigned M = stateNodes[N].size();

		for (unsigned i = 1; i < M; i++) {
			for (unsigned j = mht::kNumberOfBackSteps; j > 0; j--) {
				// Get current sepsets and scope
				emdw::RVIds pastVars = stateNodes[N-j][i]->getSepset(stateNodes[N - (j+1)][i]);
				emdw::RVIds presentVars = stateNodes[N-j][i]->getSepset( stateNodes[N - (j-1)][i] );

				// Recalibrate - marginalize onto past variables and reconstruct the joint distribution
				rcptr<Factor> pastMarginal = stateNodes[N-j][i]->marginalize(pastVars, true);
				rcptr<Factor> stateJoint = uniqptr<Factor>(new CGM( pastMarginal, 
						mht::kMotionModel, 
						presentVars,
						mht::kRCovMat  ));

				// Prune and merge
				//std::dynamic_pointer_cast<CGM>(stateJoint)->pruneAndMerge();
				stateNodes[N-j][i]->setFactor(stateJoint);

				// Determine the outgoing message using BUP
				rcptr<Factor> receivedMessage = stateNodes[N-j][i]->getReceivedMessage( stateNodes[N-(j-1)][i] );
				rcptr<Factor> outgoingMessage = 
					stateNodes[N-j][i]->marginalize(presentVars, true)->cancel(receivedMessage);

				// Receiving node absorbs and logs the message
				stateNodes[N-(j-1)][i]->inplaceAbsorb( outgoingMessage.get()  );
				stateNodes[N-(j-1)][i]->inplaceCancel( receivedMessage.get() );
				stateNodes[N-(j-1)][i]->logMessage( stateNodes[N-j][i], outgoingMessage  );
			} // for	
		} // for
	} // if
} // forwardPass()

double calculateEvidence(const unsigned N) {
	unsigned M = stateNodes[N].size();
	double odds = 0;

	// Determine the log-odds - excluding vacuous sponge.
	for (unsigned i = 0; i < M; i++) {
		std::vector<double> weights = std::dynamic_pointer_cast<CGM>(stateNodes[N][i]->getFactor())->getWeights();

		/*
		std::cout << "Target " << i << std::endl;
		std::cout << "Weights: " << weights << std::endl;
		*/

		double mass = 0;
		for (unsigned j = 0; j < weights.size(); j++) mass += weights[j];
		odds += log(mass);
	} // for

	return odds;
} // calculateEvidence()

void extractStates(const unsigned N) {
	//unsigned M = stateNodes[N].size();
	
	for (unsigned i = 1; i < 2; i++) {
		emdw::RVIds sepset = stateNodes[N][i]->getSepset( stateNodes[N - 1][i] );
		rcptr<Factor> marginal = stateNodes[N][i]->marginalize(sepset, true);

		std::vector<rcptr<Factor>> comps = std::dynamic_pointer_cast<CGM>( marginal )->getComponents();
		for (unsigned i = 0; i < comps.size(); i++) {
			ColVector<double> mean =  std::dynamic_pointer_cast<GC>(comps[i])->getMean();
			std::cout << N-1 << ";" << mean[0] << ";" << mean[2] << ";" << mean[4] << std::endl;
		} // for
	} // for	
} // extractStates()
