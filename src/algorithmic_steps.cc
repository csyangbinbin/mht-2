/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs the required tracking algorithms.
 *************************************************************************/
#include "algorithmic_steps.hpp"

void predictStates() {
	// Time step and current number of targets
	unsigned N = stateNodes.size();
	unsigned M = currentStates[N-1].size();

	// Resize for indices access
	stateNodes[N].resize(M);
	currentStates[N].resize(M);
	virtualMeasurementVars.resize(M);
	predMarginals.resize(M);

	for (unsigned i = 0; i < M; i++) {
		// Get the marginal over previous variables
		currentStates[N][i] = addVariables(variables, vecX, elementsOfX, mht::kStateSpaceDim);
		rcptr<Factor> prevMarginal = (stateNodes[N-1][i])->marginalize(elementsOfX[currentStates[N-1][i]]);

		// Create a new factor over current variables
		rcptr<Factor> stateJoint = uniqptr<Factor>(new CGM( prevMarginal, 
					initialiseMotionModel(), 
					elementsOfX[currentStates[N][i]],
					initialiseRCovMat()  ));
		stateNodes[N][i] = uniqptr<Node> (new Node(stateJoint));

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

void measurementUpdate() {
	// Time step and current number of targets
	unsigned N = stateNodes.size() - 1;
	unsigned M = currentStates[N].size();

	std::cout << N << std::endl;

	currentMeasurements[N].clear();

	// For each sensor
	for (unsigned i = 0; i < mht::kNumSensors; i++) {
		std::cout << "Sensor " << i << std::endl;

		std::vector<ColVector<double>> measurements = measurementManager->getSensorPoints(i, N + 6); // Hack time offset.

		std::vector<emdw::RVIdType> sensorMeasurements; // Local scope

		std::map<emdw::RVIdType, rcptr<DASS>> assocHypotheses; assocHypotheses.clear();
		std::map<emdw::RVIdType, ColVector<double>> colMeasurements; colMeasurements.clear();
		std::map<emdw::RVIdType, emdw::RVVals> valMeasurements; valMeasurements.clear();

		// For each measurement
		for (unsigned j = 0; j < measurements.size(); j++) { 
			if (measurements[j].size()) {
				// Assign new identity to the measurement
				emdw::RVIdType z = addVariables(variables, vecZ, elementsOfZ, mht::kMeasSpaceDim);
				sensorMeasurements.push_back(z);
				currentMeasurements[N].push_back(z);

				assocHypotheses[z] = uniqptr<DASS>(new DASS{0});
				colMeasurements[z] = measurements[j];
				valMeasurements[z] = emdw::RVVals{measurements[j][0], measurements[j][1]};

				// Form hypotheses over each measurement
				for (unsigned k = 0; k < M; k++) {
					double distance = 
						(std::dynamic_pointer_cast<GC>(validationRegion[k][j]))->mahanalobis(measurements[j]);
						
					if (distance < mht::kValidationThreshold) {
						assocHypotheses[z]->push_back(k+1); // Note the addition!
					} // if
				} // for
			} // if
		} // for

		if (assocHypotheses.size()) {
			std::map<emdw::RVIdType, rcptr<Factor>> distributions = graphBuilder->getMarginals(assocHypotheses);
			
			for (emdw::RVIdType v : sensorMeasurements) {
				DASS domain = *assocHypotheses[v];

				// Clutter is a big flat Gaussian for now.
				std::vector<rcptr<Factor>> conditionalList(domain.size());
				conditionalList[0] = uniqptr<Factor>(new CGM( elementsOfZ[v], {1.0}, {colMeasurements[v]}, {mht::kClutterCov} ));

				for (unsigned j = 1; j < domain.size(); j++) {
					unsigned p = (unsigned) (domain[j] - 1); // Note the subtraction
					
					// Create a new scope
					emdw::RVIds newScope = predMarginals[p]->getVars();
					newScope.insert(newScope.end(), elementsOfZ[v].begin(), elementsOfZ[v].end());

					// Product of predicted marginals
					conditionalList[j] = uniqptr<Factor>( predMeasurements[p][i]->copy(newScope, false) );
					conditionalList[0] = conditionalList[0]->absorb( predMarginals[p] );

					for (unsigned k = 0; k < M; k++) {
						if (k != p ) {
							conditionalList[j] = conditionalList[j]->absorb(predMarginals[k]);
						} // if
					} // for
				} // for

				rcptr<Factor> clg = uniqptr<Factor>(new LinearGaussian(distributions[v], conditionalList));

			} // for
		} // if
	} // for

	virtualMeasurementVars.clear();
	validationRegion.clear();
} // measurementUpdate()
