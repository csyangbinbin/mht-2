/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for the tentative graph builder class. The GraphBuilder
 * class creates a pairwise network of measures over the 
 * association hypotheses.
 *************************************************************************/

#ifndef GRAPH_BUILDER_HPP
#define GRAPH_BUILDER_HPP

#include <vector>
#include <random>
#include "emdw.hpp"
#include "factor.hpp"
#include "discretetable.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "lbu_cg.hpp"

/**
 * The GraphBuilder class creates a pairwise network of 
 * measures over the association hypotheses and 
 * determines the marginal beliefs held over them.
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
		 * Default destructor
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
		 * Return the marginal beliefs over the association variables.
		 */
		std::map<emdw::RVIdType, rcptr<Factor>> getMarginals(std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const;

	private:
		/**
		 * Get the association RV IDs from the given map. 
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
		 * Construct DiscreteTable factors over single association hypotheses.
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
		 * Constructs the pairwise factors required in the network, creates a
		 * cluster graph, passes messages and extracts the marginals.
		 *
		 * @param vars The association variables contianed with the map.
		 *
		 * @param assocHypotheses The association hypotheses formed over
		 * each measurement, presented as a DiscreteTables domain.
		 *
		 * @param dist A map of the association variable to the belief
		 * held over it.
		 *
		 * @return A map of the association variables to the marginal 
		 * beliefs held over them.  
		 */
		std::map<emdw::RVIdType, rcptr<Factor>> constructClusters(
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
