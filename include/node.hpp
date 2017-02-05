/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for a general cluster node.
 *************************************************************************/
#ifndef NODE_HPP
#define NODE_HPP

#include "factor.hpp"
#include "factoroperator.hpp"
#include "emdw.hpp"
#include "anytype.hpp"
#include "gausscanonical.hpp"
#include "v2vtransform.hpp"

/**
 * @brief A class representing a general cluster node.
 *
 * A class reprensenting a node in a cluster graph. A
 * wrapper type for a Factor.
 *
 * @author SCJ Robertson
 * @since 05/02/17
 */
class Node {

	public:
		/**
		 * Default constructor.
		 */
		Node();

		/**
		 * Default destructor
		 */
		~Node();

	public:
		/** 
		 * @brief Inplace normalization.
		 *
		 * Inplace normalization.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 */
		void inplaceNormalize (FactorOperator *procPtr = 0);

		/**
		 * @brief Normalization.
		 *
		 * Normalization.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 *
		 * @return A unique pointer to a normalized Factor.
		 */
		uniqptr<Factor> normalize (FactorOperator *procPtr = 0) const;

		/**
		 * @brief Inplace multiplication.
		 *
		 * Inplace multiplication. The multiplier is allowed to be a CanonicalGaussianMixture
		 * or a GaussCanonical.
		 *
		 * @param rhsPtr The divisor. Can be CanonicalGaussianMixture or
		 * GaussCanonical.
		 */
		void inplaceAbsorb (FactorOperator *procPtr = 0);

		/**
		 * @brief Multiplication.
		 *
		 * Multiplication. The multiplier is allowed to be a CanonicalGaussianMixture
		 * or a GaussCanonical.
		 *
		 * @param rhsPtr The divisor. Can be CanonicalGaussianMixture or
		 * GaussCanonical.
		 *
		 * @return A uniqptr to the product Factor.
		 */
		uniqptr<Factor> absorb (FactorOperator *procPtr = 0) const;

		/**
		 * @brief Inplace division.
		 *
		 * Inplace division by another factor.
		 *
		 * @param rhsPtr The divisor.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 */
		void inplaceCancel (FactorOperator *procPtr = 0);

		/**
		 * @brief Division.
		 *
		 * Division by another factor.
		 *
		 * @param rhsPtr The divisor.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 *
		 * @return A unique pointer to the quotient Factor.
		 */
		uniqptr<Factor> cancel (FactorOperator *procPtr = 0) const;

		/**
		 * @brief Marginalization.
		 *
		 * Marginalize out the given variables.
		 *
		 * @param variablesToKeep The variables which will not be marginalized out.
		 *
		 * @param presorted Is variablesToKeep sorted already?
		 *
		 * @param procPtr A unique pointer to some user defined FactorOperator
		 * process.
		 *
		 * @return A unique pointer to the scoped reduced Factor.
		 */
		uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep, 
				bool presorted = false, FactorOperator* procPtr = 0) const;

		/**
		I * @brief Observe and Reduce
		 *
		 * Introduce evidence and collapse the factor.
		 *
		 * @param variables The variables which have been opbsevered in some 
		 * given state.
		 *
		 * @param assignedVals The values (state) of the given variables.
		 *
		 * @param presorted Are the given variables already sorted?
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 *
		 * @return A unique pointer to resulting Factor.
		 */
		uniqptr<Factor> observeAndReduce( const emdw::RVIds& variables,
				const emdw::RVVals& assignedVals, bool presorted = false,
				FactorOperator* procPtr = 0) const;

	public:
		/**
		 * @brief Return variables.
		 */
		emdw::RVIds getVars() const;

		/**
		 * @brief Return factor.
		 */
		rcptr<Factor> getFactor() const;

		/**
		 * @brief Return the adjacent nodes
		 */
		std::vector<rcptr<Node>> getAdjacentNodes() const;

	private:
		// Current information and scope
		emdw::RVIds vars_;
		rcptr<Factor> factor_;
		
		// Neighbouring vertices
		std::vector<rcptr<Node>> adjacencyList_;
		std::map<rcptr<Node>, emdw::RVIds> sepsets_;
		
		// Past and passed information
		rcptr<Factor> prevFactor_;
		std::map<rcptr<Node>, rpctr<Factor> sentMsg_;
}; // Node

#endif // NODE_HPP
