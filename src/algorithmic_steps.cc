/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs the required tracking algorithms.
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

	for (unsigned i = 1; i < M; i++) {
		//std::cout << "Target " << stateNodes[N-1][i]->getIdentity() << std::endl;

		// Get the marginal over previous variables
		currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
		rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);
		
		// Prune and merge
		
		// Create a new factor over current variables
		rcptr<Factor> stateJoint = uniqptr<Factor>(new CGM( prevMarginal, 
					mht::kMotionModel, 
					elementsOfX[currentStates[N][i]],
					mht::kRCovMat  ));
		stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint, stateNodes[N-1][i]->getIdentity() ) );

		// Link Node to preceding node
		stateNodes[N-1][i]->addEdge(stateNodes[N][i], elementsOfX[currentStates[N-1][i]]);
		stateNodes[N][i]->addEdge(stateNodes[N-1][i], elementsOfX[currentStates[N-1][i]], prevMarginal);

		// Determine the marginal
		predMarginals[i] = stateJoint->marginalize(elementsOfX[currentStates[N][i]]);
		predMarginals[i]->inplaceNormalize();

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
	// Time step and current number of targets
	unsigned M = currentStates[N].size();

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

		// If 
		if (assocHypotheses.size()) {
			// Determine the 'prior' over association hypotheses
			std::map<emdw::RVIdType, rcptr<Factor>> distributions = graphBuilder->getMarginals(assocHypotheses);
	
			for (unsigned j = 0; j < sensorMeasurements.size(); j++) {
				emdw::RVIdType a = associations[j];
				emdw::RVIdType z = sensorMeasurements[j];

				DASS domain = *assocHypotheses[a];
				unsigned domSize = domain.size();

				//std::cout << "domains: " << domain << std::endl;
				if (domSize > 1) {
					// Create a conditional map for the ConditionalGaussian
					std::map<emdw::RVIdType, rcptr<Factor>> conditionalList; conditionalList.clear();

					for (unsigned k = 1; k < domSize; k++) {
						unsigned p = domain[k]; 
						
						// Create a new scope
						emdw::RVIds newScope = predMarginals[p]->getVars();
						newScope.insert(newScope.end(), elementsOfZ[z].begin(), elementsOfZ[z].end());

						// Product of predicted marginals
						conditionalList[p] = uniqptr<Factor>( predMeasurements[p][i]->copy(newScope, false) );

						// Clutter needs to be a joint over all 'included' predicted states
						if (k == 1) conditionalList[0] = uniqptr<Factor>(predMarginals[p]->copy());
						else conditionalList[0] = conditionalList[0]->absorb(predMarginals[p]);

						// Multiply the predicted states by the likelihood function
						for (unsigned l = 1; l < domSize; l++) {
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
			
					// Adjust for a uniform clutter rate over sensor volume
					std::dynamic_pointer_cast<CGM>(conditionalList[0])->adjustMass(1e-3);

					// Create ConditionalGauss - observeAndReduce does work, but for scope reasons this is easier.
					rcptr<Factor> clg = uniqptr<Factor>(new CLG(distributions[a], conditionalList));

					// Create a measurement node and connect it to state nodes
					rcptr<Node> measNode = uniqptr<Node>(new Node(clg));
					for (unsigned k = 1; k < domain.size(); k++) {
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
	//std::cout << "Measurement updates M: " << M << std::endl;

	for (unsigned i = 0; i < M; i++) {
		std::vector<std::weak_ptr<Node>> adjacent = measurementNodes[N][i]->getAdjacentNodes();

		for (unsigned j = 0; j < adjacent.size(); j++) {
			// Get the neighbouring state node and message it sent to the measurement clique
			rcptr<Node> stateNode = adjacent[j].lock();
			rcptr<Factor> receivedMessage = measurementNodes[N][i]->getReceivedMessage( stateNode );

			// Determine the outgoing message
			emdw::RVIds sepset = measurementNodes[N][i]->getSepset( stateNode );
			rcptr<Factor> outgoingMessage =  measurementNodes[N][i]->marginalize( sepset, true )->cancel(receivedMessage);

			// Get the factor, absorb  the incoming message and then prune
			rcptr<Factor> factor = stateNode->getFactor();
			factor->inplaceAbsorb( outgoingMessage );
			std::dynamic_pointer_cast<CGM>(factor)->pruneAndMerge();

			// Update the factor
			stateNode->setFactor(factor);
			stateNode->logMessage( measurementNodes[N][i] ,outgoingMessage );
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
				rcptr<Factor> outgoingMessage = stateNodes[N-j][i]->marginalize(sepset, true)->cancel(receivedMessage);

				// Received node absorbs and logs the message
				stateNodes[N-(j+1)][i]->absorb( outgoingMessage.get()  );
				stateNodes[N-(j+1)][i]->logMessage( stateNodes[N-j][i], outgoingMessage  );
			} // for

			// Forward pass
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
				std::dynamic_pointer_cast<CGM>(stateJoint)->pruneAndMerge();
				stateNodes[N-j][i]->setFactor(stateJoint);

				// Determine the outgoing message using BUP
				rcptr<Factor> receivedMessage = stateNodes[N-j][i]->getReceivedMessage( stateNodes[N-(j-1)][i] );
				rcptr<Factor> outgoingMessage = 
					stateNodes[N-j][i]->marginalize(presentVars, true)->cancel(receivedMessage);

				// Receiving node absorbs and logs the message
				stateNodes[N-(j-1)][i]->absorb( outgoingMessage.get()  );
				stateNodes[N-(j-1)][i]->logMessage( stateNodes[N-j][i], outgoingMessage  );
			} // for

			// TODO: Remove this and make a logger at some point
			if ( stateNodes[N][i]->getIdentity() == 1 ) {
				emdw::RVIds sepset = stateNodes[N][i]->getSepset( stateNodes[N - 1][i] );
				rcptr<Factor> marginal = stateNodes[N][i]->marginalize(sepset, true);

				std::vector<rcptr<Factor>> comps = std::dynamic_pointer_cast<CGM>( marginal )->getComponents();
				for (unsigned i = 0; i < comps.size(); i++) {
					ColVector<double> mean =  std::dynamic_pointer_cast<GC>(comps[i])->getMean();
					std::cout << N-1 << ";" << mean[0] << ";" << mean[2] << ";" << mean[4] << std::endl;
				} // for
			} // if


		} // for
	} // if

	stateNodes[N-mht::kNumberOfBackSteps].clear();
	measurementNodes[N-mht::kNumberOfBackSteps].clear();
} // smoothTrajectory()

void modelSelection(unsigned const N) {

} // modelSelection()
