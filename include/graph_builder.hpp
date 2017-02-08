/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for the tentative graph builder class. The GraphBuilder
 * class creates a pairwise network of unnormalized measures over the 
 * association hypotheses.
 *************************************************************************/

#ifndef GRAPH_BUILDER_HPP
#define GRAPH_BUILDER_HPP

#include <vector>
#include "emdw.hpp"
#include "factor.hpp"
#include "discretetable.hpp"
#include "graph.hpp"

// Forward declaration
class GraphBuilder;

/**
 * The GraphBuilder class creates a pairwise network of unnormalized 
 * measures over the association hypotheses. 
 *
 * This does not construct a general cluster graph, but the name
 * stuck.
 *
 * @author SCJ Robertson
 * @since 23/01/17
 */
class GraphBuilder {

	private:
		typedef unsigned short T;
		typedef DiscreteTable<T> DT;
		typedef std::vector<T> DASS;

	public:
		/**
		 * Construct a pairwise network of unnormalized measures over
		 * the association hypotheses.
		 *
		 * @param assoc_hypotheses The association hypotheses formed
		 * over each measurement. Presented as a DiscreteTable's domain.
		 */
		GraphBuilder(
			const std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses,
			const double floor = 0.0,
			const double margin = 0.0,
			const double defProb = 0.0,
			const rcptr<FactorOperator>& inplaceNormalizer = 0,
			const rcptr<FactorOperator>& normalizer = 0,
			const rcptr<FactorOperator>& marginalizer = 0
			);

		/**
		 * Default constructor.
		 */
		GraphBuilder();

		/**
		 * Default destructor.
		 */
		~GraphBuilder();

	private:
		/**
		 * Get the association RV IDs from the given map. Initialises
		 * the a_ member and sorts the association hypotheses domains.
		 */
		void extractRVIds();

		/**
		 * Construct DiscreteTable factors over single
		 * association hypotheses. Initialises the distribution_
		 * member.
		 */
		void constructDistributions();

		/**
		 * Constructs the pairwise factors required in the network.
		 * Initialises the cluster_ member.
		 *
		 * TODO: Better implementation? Pairwise connections are
		 * brute force.
		 */
		void constructClusters ();


	private:
		emdw::RVIds a_;
		std::vector<rcptr<Graph>> graphs_;

		std::map <emdw::RVIdType, rcptr<Factor>> dist_;
		std::map <emdw::RVIdType, rcptr<DASS>> assocHypotheses_;

		rcptr<FactorOperator> inplaceNormalizer_;
		rcptr<FactorOperator> normalizer_;
		rcptr<FactorOperator> marginalizer_;

		double floor_;
		double margin_;
		double defProb_;
}; // GraphBuilder()

#endif // GRAPH_BUILDER_HPP
