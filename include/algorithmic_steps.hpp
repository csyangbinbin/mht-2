/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file declaring function prototypes for the algorithmic
 * steps used in main.cc
 *************************************************************************/
#ifndef ALGORITHMICSTEPS_HPP
#define ALGORITHMICSTEPS_HPP

#include <iostream>
#include <map>
#include "emdw.hpp"
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "conditional_gaussian.hpp"
#include "transforms.hpp"
#include "utils.hpp"
#include "system_constants.hpp"

/**
 * @brief Predict the current state of the x variables
 *
 * @param N The current time index.
 */
void predictStatesSU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion);

/**
 * @brief Predict the current state of the x variables
 *
 * @param N The current time index.
 */
void predictStatesAU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes);

/**
 * @brief Forms hypotheses and creates current measurement distributions.
 *
 * @param N The current time index.
 */
void createMeasurementDistributionsSU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion);

/**
 * @brief Forms hypotheses and creates current measurement distributions.
 *
 * @param N The current time index.
 */
void createMeasurementDistributionsAU(const unsigned N,
		const unsigned sensorNumber,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion);

/**
 * @brief Performs measurement update on exisitng targets.
 *
 * @param N The current time index.
 */
void measurementUpdateSU(const unsigned N,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes);

/**
 * @brief Forms hypotheses and creates current measurement distributions.
 *
 * @param N The current time index.
 */
void measurementUpdateAU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion);


/**
 * @brief Smoothes existing targets trajectories.
 *
 * @param N The current time index.
 */
void smoothTrajectory(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes);

/**
 * @brief Decide whether to add new targets.
 *
 * @param N The current time index.
 */
void modelSelectionSU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion
		);

/**
 * @brief Decide whether to add new targets.
 *
 * @param N The current time index.
 */
void modelSelectionAU(const unsigned N,
		std::map<unsigned, emdw::RVIds>& currentStates,
		emdw::RVIds& virtualMeasurementVars, 
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes,
		std::map<unsigned, std::vector<rcptr<Node>>>& measurementNodes,
		std::vector<rcptr<Factor>>& predMarginals,
		std::map<unsigned, std::vector<rcptr<Factor>>>& predMeasurements,
		std::map<unsigned, std::vector<rcptr<Factor>>>& validationRegion
		);

/**
 * @brief Propagates new targets forward, possibly
 * reintroducing new evidence.
 *
 * @param N The current time index.
 */
void forwardPass(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes);

/**
 * @brief Remove all targets which have grounded.
 *
 * @param N The current time index
 */
void removeStates(const unsigned N, 
		std::map<unsigned, emdw::RVIds>& currentStates,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes);

#endif // ALGORITHMICSTEPS_HPP
