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
#include <random>
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
 * TODO: This should be a module that just shits out
 * association Graphs for given hypotheses. I need to
 * rewrite most of this object.
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
		 */
		GraphBuilder(
			const double floor = 0.0,
			const double margin = 0.0,
			const double defProb = 0.0,
			const rcptr<FactorOperator>& inplaceNormalizer = 0,
			const rcptr<FactorOperator>& normalizer = 0,
			const rcptr<FactorOperator>& marginalizer = 0
			);

		/**
		 * Default destructor.
		 */
		~GraphBuilder();

	public:
		/**
		 * Return the constructed graphs.
		 */
		std::vector<rcptr<Graph>> getGraphs(std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const;

	private:
		/**
		 * Get the association RV IDs from the given map. Initialises
		 * the a_ member and sorts the association hypotheses domains.
		 *
		 * @param assocHypotheses The association hypotheses formed over
		 * each measurement, presented as a DiscreteTables domain.
		 *
		 * @return emdw::RVIds of the variables contained in the map. 
		 */
		emdw::RVIds extractRVIds(
				const std::map<emdw::RVIdType, 
				rcptr<DASS>>& assocHypotheses) const;

		/**
		 * Construct DiscreteTable factors over single
		 * association hypotheses. Initialises the distribution_
		 * member.
		 *
		 * @param vars The association variables contianed with the map.
		 *
		 * @param assocHypotheses The association hypotheses formed over
		 * each measurement, presented as a DiscreteTables domain.
		 *
		 * @return A map of the association variable to the distribution 
		 * held over it.
		 */
		std::map<emdw::RVIdType, rcptr<Factor>> constructDistributions(
				const emdw::RVIds& vars,
				std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const;

		/**
		 * Constructs the pairwise factors required in the network.
		 * Initialises the cluster_ member.
		 *
		 * TODO: This is the worst code ever. May take longer to construct
		 * than explicitly using a massive table.
		 *
		 * @return A vector of disjoint graphs, representing 
		 * association networks.
		 */
		std::vector<rcptr<Graph>> constructClusters(
				const emdw::RVIds& vars,
				std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses,
				std::map<emdw::RVIdType, rcptr<Factor>>& dist
				) const;

	private:
		// DiscreteTable properties
		double floor_;
		double margin_;
		double defProb_;

		// DiscreteTable factor operators
		rcptr<FactorOperator> inplaceNormalizer_;
		rcptr<FactorOperator> normalizer_;
		rcptr<FactorOperator> marginalizer_;

}; // GraphBuilder()

#endif // GRAPH_BUILDER_HPP
