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
		GraphBuilder(const std::map<emdw::RVIdType, rcptr<DASS>>& assoc_hypotheses);

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
		 * association hypotheses.
		 */
		void constructDistributions();

		/**
		 * Constructs the pairwise factors required in the network.
		 */
		void constructClusters ();

		/**
		 * Returns the RVIds of the association hypotheses.
		 *
		 * @return The IDs of the association hypotheses.
		 */
		emdw::RVIds getRVIds() const;

	private:
		emdw::RVIds rv_ids_;
		std::vector <rcptr<Factor>> distributions_;
		std::vector <rcptr<Factor>> cluster_;


		std::map <emdw::RVIdType, rcptr<DASS>> assoc_hypotheses_;
		std::map <Idx2, emdw::RVIds> sepvecs_;
};

#endif // GRAPH_BUILDER.HPP
