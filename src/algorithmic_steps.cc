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

		// Moment matching - Horrible, but nothing else seems to work.
		rcptr<Factor> momentMatched = std::dynamic_pointer_cast<CGM>(prevMarginal)->momentMatch();
		prevMarginal = uniqptr<Factor> (new CGM(momentMatched->getVars(), {momentMatched} ));

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

		if ( stateNodes[N-1][i]->getIdentity() == 1 ) {
			std::vector<rcptr<Factor>> comps = std::dynamic_pointer_cast<CGM>( predMarginals[i] )->getComponents();
			for (rcptr<Factor> c: comps) {
				ColVector<double> mean =  std::dynamic_pointer_cast<GC>(c)->getMean();
				std::cout << N << ";" << mean[0] << ";" << mean[2] << ";" << mean[4] << std::endl;
			}
		}

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

	// Clear
	stateNodes[N-1].clear();
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
					// Clutter distribution is a big flat Gaussian for now.
					std::map<emdw::RVIdType, rcptr<Factor>> conditionalList; conditionalList.clear();
					conditionalList[0] = uniqptr<Factor>(new CGM( 
								elementsOfZ[z], 
								{1.0}, 
								{colMeasurements[z]}, 
								{mht::kClutterCov} ));

					for (unsigned k = 1; k < domSize; k++) {
						unsigned p = domain[k]; 
						
						// Create a new scope
						emdw::RVIds newScope = predMarginals[p]->getVars();
						newScope.insert(newScope.end(), elementsOfZ[z].begin(), elementsOfZ[z].end());

						// Product of predicted marginals
						conditionalList[p] = uniqptr<Factor>( predMeasurements[p][i]->copy(newScope, false) );
						conditionalList[0] = conditionalList[0]->absorb( predMarginals[p] );

						for (unsigned l = 1; l < domSize; l++) {
							unsigned q = domain[l];
							if ( q != p ) conditionalList[p] = 
								conditionalList[p]->absorb(predMarginals[q] );
						} // for
						conditionalList[p]->inplaceNormalize();
					} // for
					conditionalList[0]->inplaceNormalize();

					// Create and reduce a CLG
					rcptr<Factor> clg = uniqptr<Factor>(new LinearGaussian(distributions[a], conditionalList));
					rcptr<Factor> reduced = clg->observeAndReduce(elementsOfZ[z], 
						emdw::RVVals{colMeasurements[z][0], colMeasurements[z][1]}, 
						true);

					// Create a measurement node and connect it to state nodes
					rcptr<Node> measNode = uniqptr<Node>(new Node(reduced));
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

	for (unsigned i = 0; i < M; i++) {
		std::vector<std::weak_ptr<Node>> adjacent = measurementNodes[N][i]->getAdjacentNodes();

		for (unsigned j = 0; j < adjacent.size(); j++) {
			// Get the neighbouring state node and message it sent to the measurement clique
			rcptr<Node> stateNode = adjacent[j].lock();
			rcptr<Factor> receivedMessage = measurementNodes[N][i]->getReceivedMessage( stateNode );

			// Determine the outgoing message
			emdw::RVIds sepset = measurementNodes[N][i]->getSepset( stateNode );
			rcptr<Factor> outgoingMessage =  measurementNodes[N][i]->marginalize( sepset, true )->cancel(receivedMessage);

			// The state node absorbs the incoming message
			stateNode->inplaceAbsorb( outgoingMessage.get() );
			stateNode->logMessage( measurementNodes[N][i] ,outgoingMessage );

		} // for
	} // for
	measurementNodes[N].clear();
} // measurementUpdate()
