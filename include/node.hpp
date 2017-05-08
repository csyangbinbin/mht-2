/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for a general cluster node.
 *************************************************************************/
#ifndef NODE_HPP
#define NODE_HPP

#include <iostream>
#include "factor.hpp"
#include "factoroperator.hpp"
#include "emdw.hpp"
#include "anytype.hpp"

// Forward Declaration
class Node;

/**
 * @brief A class representing a general cluster node.
 *
 * A class reprensenting a node in a cluster graph. A
 * wrapper type for a Factor, this uses an
 * adjacency list representation for a Graph.
 *
 * @author SCJ Robertson
 * @since 05/02/17
 */
class Node {

	public:
		/**
		 * @brief Node constructor.
		 *
		 * Create a general cluster node given the factor it
		 * contains and its adjacent nodes.
		 * 
		 * @param factor The cluster the node is to contain.
		 */
		Node(const rcptr<Factor>& factor, const unsigned N = 0);

		/**
		 * @brief Default destructor.
		 */
		~Node();

	public:
		/**
		 * @brief Add an edge.
		 *
		 * Connect this cluster node to an
		 * adjacent cluster node.
		 *
		 * @param w An adjacent cluster node.
		 *
		 * @param sepset The sepset variables shared between
		 * the two nodes.
		 *
		 * @param message An initial message sent to the cluster.
		 */
		void addEdge(const rcptr<Node>& w, const emdw::RVIds& sepset, const rcptr<Factor>& message = 0);

		/**
		 * @brief Remove and edge.
		 *
		 * Remove an edge between two clusters.
		 *
		 * @param w An adjacent cluster node.
		 */
		void removeEdge(const rcptr<Node>& w);

		/**
		 * @brief Add a newly received message to the log.
		 *
		 * Add a newly received message to the message 
		 * log.
		 *
		 * @param w The neighbouring node in the graph that
		 * sent the message.
		 *
		 * @param message The newly received message.
		 */
		void logMessage(const rcptr<Node>& w, const rcptr<Factor>& message);

		/**
		 *
		 */
		void setFactor(const rcptr<Factor>& factor);

		/**
		 * @brief Cache a the current status of the factor.
		 * 
		 * Cache the current status of the factor. This allows for
		 * convergence tests during the next round of message passing.
		 *
		 * @param factor The current status of the factor.
		 */
		void cacheFactor(const rcptr<Factor>& factor);

	public:
		/**
		 * @brief Return a unique identity.
		 */
		unsigned getIdentity() const;

		/**
		 * @brief Return variables.
		 */
		emdw::RVIds getVars() const;

		/**
		 * @brief Return factor.
		 */
		rcptr<Factor> getFactor() const;

		/**
		 * @brief Return the previous the cached factor
		 */
		rcptr<Factor> getCachedFactor() const;

		/**
		 * @brief Return sepset variables
		 */
		emdw::RVIds getSepset(const rcptr<Node>& w);

		/**
		 * @brief Return last message received from a given neighbour.
		 */
		rcptr<Factor> getReceivedMessage(const rcptr<Node>& w);

		/**
		 * @brief Return the adjacent nodes
		 */
		std::vector<std::weak_ptr<Node>> getAdjacentNodes() const;
	
	public:
		/** 
		 * @brief Inplace normalization.
		 *
		 * Inplace normalization.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 */
		void inplaceNormalize (FactorOperator* procPtr = 0);

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
		uniqptr<Factor> normalize (FactorOperator* procPtr = 0) const;

		/**
		 * @brief Inplace multiplication.
		 *
		 * Inplace multiplication. The multiplier is allowed to be a CanonicalGaussianMixture
		 * or a GaussCanonical.
		 *
		 * @param rhsPtr The divisor. Can be CanonicalGaussianMixture or
		 * GaussCanonical.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 */
		void inplaceAbsorb (const Factor* rhsPtr, FactorOperator* procPtr = 0);

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
		uniqptr<Factor> absorb (const Factor* rhsPtr, FactorOperator* procPtr = 0) const;

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
		void inplaceCancel (const Factor* rhsPtr, FactorOperator *procPtr = 0);

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
		uniqptr<Factor> cancel (const Factor* rhsPtr, FactorOperator* procPtr = 0) const;

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
		 * @brief Inplace Observe and Reduce
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
		 */
		void inplaceObserveAndReduce( const emdw::RVIds& variables,
				const emdw::RVVals& assignedVals, bool presorted = false,
				FactorOperator* procPtr = 0);

		/**
		 * @brief Observe and Reduce
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
		 * @brief Write information to an output stream.
		 */
		friend std::ostream& operator<<(std::ostream& file, const Node& node);

	private:
		// Current information and scope
		rcptr<Factor> factor_;
		unsigned N_;
		emdw::RVIds vars_;
		
		// Neighbouring vertices
		mutable std::map<std::weak_ptr<Node>, emdw::RVIds, std::owner_less<std::weak_ptr<Node>>> sepsets_;
		mutable std::vector<std::weak_ptr<Node>> adjacent_;
		
		// Past and passed information
		rcptr<Factor> prevFactor_;
		mutable std::map<std::weak_ptr<Node>, rcptr<Factor>, std::owner_less<std::weak_ptr<Node>>> recMsg_;
}; // Node

#endif // NODE_HPP
