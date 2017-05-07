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
	std::cout << "\npredictStates(" << N << ")" << std::endl;
	// Time step and current number of targets
	unsigned M = stateNodes[N-1].size();

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
	std::cout << "\n\n" << std::endl;

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
	std::cout << "\ncreateMeasurementDistributions()" << std::endl;
	
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

	std::cout << "\n\n" << std::endl;

	// Clear temporary values
	predMarginals.clear();
	virtualMeasurementVars.clear();
	validationRegion.clear();
} // createMeasurementDistributions()

void measurementUpdate(const unsigned N,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes
		) {
	unsigned M = measurementNodes[N].size();

	std::cout << "\nmeasurementUpdate(), M: " << M << std::endl;

	for (unsigned i = 0; i < M; i++) {
		std::vector<std::weak_ptr<Node>> adjacent = measurementNodes[N][i]->getAdjacentNodes();

		for (unsigned j = 0; j < adjacent.size(); j++) {
			// Get the neighbouring state node and message it sent to the measurement clique
			rcptr<Node> stateNode = adjacent[j].lock();
			std::cout << "Target: " << stateNode->getIdentity() << std::endl;

			rcptr<Factor> receivedMessage = measurementNodes[N][i]->getReceivedMessage( stateNode );
			
			// Determine the outgoing message
			emdw::RVIds sepset = measurementNodes[N][i]->getSepset( stateNode );
			rcptr<Factor> outgoingMessage =  (measurementNodes[N][i]->marginalize( sepset, true));

			// Get the factor, absorb  the incoming message and then prune
			rcptr<Factor> factor = stateNode->getFactor();
			factor->inplaceAbsorb( outgoingMessage );
			factor->inplaceCancel( receivedMessage );

			// Update the factor
			stateNode->setFactor(factor);
			stateNode->logMessage(measurementNodes[N][i], outgoingMessage);
		} // for
	} // for
	std::cout << "\n\n" << std::endl;
} // measurementUpdate()

void smoothTrajectory(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	if (N > mht::kNumberOfBackSteps) { 

		std::cout << "\nsmoothTrajectory()" << std::endl;

		unsigned M = stateNodes[N].size();

		for (unsigned i = 1; i < M; i++) {
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
		std::cout << "\n\n" << std::endl;
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
	if ( N > mht::kNumberOfBackSteps ) {

		std::cout << "modelSelection()" << std::endl;
		double modelOneOdds = calculateEvidence(N, stateNodes);

		std::cout << "modelOneOdds: " << exp(modelOneOdds) << std::endl;

		std::cout << "\n\n" << std::endl;
	} // if
} // modelSelection()

void forwardPass(unsigned const N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	if (N > mht::kNumberOfBackSteps) {
		unsigned M = stateNodes[N].size();

		std::cout << "\nforwardPass()" << std::endl;

		for (unsigned i = 1; i < M; i++) {
			for (unsigned j = mht::kNumberOfBackSteps; j > 0; j--) {
				std::cout << "j: " << j << ", N-(j-1): " << N-(j-1) << std::endl;

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
		std::cout << "\n\n" << std::endl;
	} // if
} // forwardPass()

double calculateEvidence(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {

	std::cout << "calculateEvidence()" << std::endl;

	unsigned M = stateNodes[N].size();
	double odds = 0;

	// Determine the log-odds - excluding vacuous sponge.
	for (unsigned i = 0; i < M; i++) {
		double mass = std::dynamic_pointer_cast<CGM>(stateNodes[N][i]->getFactor())->getMass();
		odds += log(mass);
	} // for

	return odds;
} // calculateEvidence()

void extractStates(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	std::cout << "\nextractStates()" << std::endl;
	unsigned M = stateNodes[N].size();
	
	for (unsigned i = 1; i < M; i++) {
		rcptr<Factor> marginal = stateNodes[N][i]->marginalize(elementsOfX[currentStates[N][i]], true);
		std::vector<rcptr<Factor>> comps = std::dynamic_pointer_cast<CGM>( marginal )->getComponents();
		
		for (unsigned j = 0; j < comps.size(); j++) {
			ColVector<double> mean =  std::dynamic_pointer_cast<GC>(comps[j])->getMean();
			double mass = std::dynamic_pointer_cast<GC>(comps[j])->getMass();
			std::cout << N << ";" << i << ";" << mean[0] << ";" << mean[2] << ";" << mean[4] 
				<< "; Mass: " << mass << std::endl;
		} // for
	} // for

	std::cout << "\n\n" << std::endl;
} // extractStates()
