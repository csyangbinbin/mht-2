/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * A library of useful helper functions which don't really fit in anywhere 
 * else. 
 *************************************************************************/
#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <map>
#include "emdw.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "node.hpp"
#include "system_constants.hpp"

/**
 * @brief Add a new vector valued variable to map collection.
 *
 * Add a new vector valued variable to map collection.
 *
 * @param globalVariables A vector of all exisitng variables.
 *
 * @param localVariables A vector containing IDs of all
 * specific types of vector valued variables.
 *
 * @param map A map of vector valued IDs to its elements
 *
 * @param L The dimension of the vector variable.
 *
 * @return The identity of the new vector valued
 * variable.
 */
unsigned addVariables (emdw::RVIds& globalVariables, 
		emdw::RVIds& localVariables, 
		std::map<unsigned,emdw::RVIds>& map,
		const unsigned L);
/**
 * @brief Determine the evidence, after
 * smoothing and measurement update.
 *
 * @param N The current time index.
 *
 * @return The evidence provided in logarithmic form.
 */
double calculateEvidence(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes);

/**
 * @brief Extract the targets' states at each time step.
 *
 * @param N The current time index.
 */
void extractStates(const unsigned N, 
		std::map<unsigned, emdw::RVIds>& currentStates,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes);

/**
 * @brief Determine the factorial of an unsigned integer.
 *
 * @param N An unsigned integer.
 *
 * @return The factorial of N.
 */
unsigned factorial (const unsigned N);

#endif // UTILS_HPP
