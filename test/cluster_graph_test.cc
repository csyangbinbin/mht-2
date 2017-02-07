/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for graph.hpp and node.hpp.
 *************************************************************************/
#include <iostream>
#include <map>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "anytype.hpp"
#include "emdw.hpp"
#include "oddsandsods.hpp"
#include "discretetable.hpp"
#include "graph.hpp"
#include "node.hpp"

using namespace emdw;
using namespace std;

class GraphTest : public testing::Test {
	protected:
		typedef unsigned short T;
		typedef DiscreteTable<T> DT;
		typedef std::vector<T> DASS;

	protected:
		virtual void SetUp() {
			margPtr_ = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);
			inormPtr_ = uniqptr<FactorOperator> (new DiscreteTable_InplaceNormalize<T>);
			normPtr_ = uniqptr<FactorOperator> (new DiscreteTable_MaxNormalize<T>);

			// Probability table for a0
			a0Dom_ = uniqptr<DASS>(new DASS{0, 1});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{1}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a0}, {a0Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, margPtr_, inormPtr_, normPtr_) ) );

			// Probability table for a1
			a1Dom_ = uniqptr<DASS>(new DASS{0, 1, 2});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{1}] = sparseProbs_[DASS{2}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a1}, {a1Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, margPtr_, inormPtr_, normPtr_) ) );

			// Probability table for a2
			a2Dom_ = uniqptr<DASS>(new DASS{0, 2});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{2}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a2}, {a2Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, margPtr_, inormPtr_, normPtr_) ) );
		}

		virtual void TearDown() {}

	protected:
		enum{a0, a1, a2};
		const double kFloor_ = 0.0;
		const double kMargin_ = 0.0;
		const double kDefProb_ = 0.0;

		rcptr<FactorOperator> margPtr_;
		rcptr<FactorOperator> inormPtr_;
		rcptr<FactorOperator> normPtr_;

		rcptr<DASS> a0Dom_;
		rcptr<DASS> a1Dom_;
		rcptr<DASS> a2Dom_;

		std::vector<rcptr<Factor>> hypotheses_;
		std::map<DASS, FProb> sparseProbs_;
};

TEST_F (GraphTest, ControlledExample) {
	// Create clusters
	std::vector<rcptr<Factor>> cluster;
	cluster.push_back(hypotheses_[0]->absorb(hypotheses_[1]));
	cluster.push_back(hypotheses_[0]->absorb(hypotheses_[2]));
	cluster.push_back(hypotheses_[1]->absorb(hypotheses_[2]));

	// Remove any shared hypotheses
	rcptr<DT> change = std::dynamic_pointer_cast<DT>(cluster[0]);
	change->setEntry(emdw::RVIds{a0, a1}, emdw::RVVals{T(0), T(0)}, 0);
	change->setEntry(emdw::RVIds{a0, a1}, emdw::RVVals{T(1), T(1)}, 0);
	cluster[0]->inplaceNormalize();

	change = std::dynamic_pointer_cast<DT>(cluster[1]);
	change->setEntry(emdw::RVIds{a0, a2}, emdw::RVVals{T(0), T(0)}, 0);
	cluster[1]->inplaceNormalize();
	
	change = std::dynamic_pointer_cast<DT>(cluster[2]);
	change->setEntry(emdw::RVIds{a1, a2}, emdw::RVVals{T(0), T(0)}, 0);
	change->setEntry(emdw::RVIds{a1, a2}, emdw::RVVals{T(2), T(2)}, 0);
	cluster[2]->inplaceNormalize();

	// Create the nodes
	rcptr<Node> v = uniqptr<Node>(new Node(cluster[0]));
	rcptr<Node> w = uniqptr<Node>(new Node(cluster[1]));
	rcptr<Node> x = uniqptr<Node>(new Node(cluster[2]));

	// Create the graph
	rcptr<Graph> graph = uniqptr<Graph>(new Graph());
	graph->addEdge(v, w);
	graph->addEdge(w, x);
	graph->addEdge(v, x);

	// Run a DFS message passing scheme
	for (unsigned i = 0; i < 100; i++) graph->depthFirstSearch();

	std::cout << "====================" << std::endl;
	std::cout << *(v->getFactor()) << std::endl;
	std::cout << "====================" << std::endl;
	std::cout << *(w->getFactor()) << std::endl;
	std::cout << "====================" << std::endl;
	std::cout << *(x->getFactor()) << std::endl;
}
