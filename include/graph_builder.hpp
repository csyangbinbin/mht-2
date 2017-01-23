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

	public:
		/**
		 * Default constructor.
		 *
		 * @param assoc_hypotheses Association hypotheses formed over each
		 * measurement. Provided as a DiscreteTable over a single variable a_i.
		 */
		GraphBuilder(const std::vector <rcptr<Factor>>& assoc_hypotheses);

	private:
		typedef unsigned short T;
		typedef DiscreteTable<T> DT;
		typedef std::vector<DT> DASS;

	private:
		std::vector <rcptr<Factor>> assoc_;
};

#endif // GRAPH_BUILDER.HPP
