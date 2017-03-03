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
	stateNodes[N].resize(M); stateNodes[N][0] = 0;
	currentStates[N].resize(M);
	virtualMeasurementVars.resize(M);
	predMarginals.resize(M);

	for (unsigned i = 1; i < M; i++) {
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

	currentMeasurements[N].clear();

	// For each sensor
	for (unsigned i = 0; i < mht::kNumSensors; i++) {
		std::cout << "Sensor " << i << std::endl;

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

				associations.push_back(a);
				sensorMeasurements.push_back(z);
				currentMeasurements[N].push_back(z);

				assocHypotheses[a] = uniqptr<DASS>(new DASS{0});
				colMeasurements[z] = measurements[j];

				// Form hypotheses over each measurement
				for (unsigned k = 1; k < M; k++) {
					double distance = 
						(std::dynamic_pointer_cast<GC>(validationRegion[k][j]))->mahanalobis(measurements[j]);
						
					if (distance < mht::kValidationThreshold) {
						assocHypotheses[a]->push_back(k); 
					} // if
				} // for
			} // if
		} // for

		if (assocHypotheses.size()) {
			std::map<emdw::RVIdType, rcptr<Factor>> distributions = graphBuilder->getMarginals(assocHypotheses);
	
			for (unsigned j = 0; j < sensorMeasurements.size(); j++) {
				emdw::RVIdType a = associations[j];
				emdw::RVIdType z = sensorMeasurements[j];

				// Clutter is a big flat Gaussian for now.
				std::map<emdw::RVIdType, rcptr<Factor>> conditionalList; conditionalList.clear();
				conditionalList[0] = 
					uniqptr<Factor>(new CGM( elementsOfZ[z], {1.0}, {colMeasurements[z]}, {mht::kClutterCov} ));

				DASS domain = *assocHypotheses[a];
				for (unsigned k = 1; k < domain.size(); k++) {
					unsigned p = (unsigned) domain[k]; 
					
					// Create a new scope
					emdw::RVIds newScope = predMarginals[p]->getVars();
					newScope.insert(newScope.end(), elementsOfZ[z].begin(), elementsOfZ[z].end());

					// Product of predicted marginals
					conditionalList[p] = uniqptr<Factor>( predMeasurements[p][i]->copy(newScope, false) );
					conditionalList[0] = conditionalList[0]->absorb( predMarginals[p] );

					for (unsigned l = 1; l < M; l++) {
						if (l != p ) {
							conditionalList[p] = conditionalList[p]->absorb(predMarginals[l]);
						} // if
					} // for
				} // for

				rcptr<Factor> clg = uniqptr<Factor>(new LinearGaussian(distributions[a], conditionalList));
				rcptr<Factor> reduced = clg->observeAndReduce(elementsOfZ[z], 
						emdw::RVVals{colMeasurements[z][0], colMeasurements[z][1]}, 
						true);
				rcptr<Factor> cgm = reduced->marginalize(std::dynamic_pointer_cast<LG>(reduced)->getContinuousVars());
				cgm->inplaceNormalize();

				std::vector<rcptr<Factor>> comps = std::dynamic_pointer_cast<CGM>(cgm)->getComponents();

				for (rcptr<Factor> c : comps) {
					std::cout << std::dynamic_pointer_cast<GC>(c)->getMass() << std::endl;
				}
			} // for
		} // if
	} // for

	virtualMeasurementVars.clear();
	validationRegion.clear();
} // measurementUpdate()
