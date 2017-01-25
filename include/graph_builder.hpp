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

/**
 * The GraphBuilder class creates a pairwise network of unnormalized 
 * measures over the association hypotheses.
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
		GraphBuilder(const std::map<emdw::RVIdType, rcptr<DASS>>& assoc_hypotheses,
			const rcptr<FactorOperator>& marg_ptr, 
			const rcptr<FactorOperator>& inorm_ptr,
			const rcptr<FactorOperator>& norm_ptr,
			const double floor,
			const double margin,
			const double def_prob);

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
		 * Construct DiscreteTable factors over single
		 * association hypotheses. Initialises the distribution_
		 * member.
		 */
		void ConstructDistributions();

		/**
		 * Constructs the pairwise factors required in the network.
		 * Initialises the cluster_ member.
		 *
		 * TODO: Better implementation? Pairwise connections are
		 * brute force.
		 */
		void ConstructClusters ();

		/**
		 * Get the association RV IDs from the given map. Initialises
		 * the a_ member and sorts the association hypotheses domains.
		 */
		void AcquireRVIds();

	private:
		emdw::RVIds a_;
		std::vector <rcptr<Factor>> cluster_;

		std::map <emdw::RVIdType, rcptr<Factor>> distribution_;
		std::map <emdw::RVIdType, rcptr<DASS>> assoc_hypotheses_;
		std::map <Idx2, emdw::RVIds> sepvecs_;

		rcptr<FactorOperator> marg_ptr_;
		rcptr<FactorOperator> inorm_ptr_;
		rcptr<FactorOperator> norm_ptr_;

		double floor_;
		double margin_;
		double def_prob_;
};

#endif // GRAPH_BUILDER.HPP
