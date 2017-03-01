/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for the loop assocation tests
 *************************************************************************/
#include <iostream>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "anytype.hpp"
#include "emdw.hpp"
#include "vecset.hpp"
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "clustergraph.hpp"
#include "lbp_cg.hpp"
#include "lbu_cg.hpp"
#include "graph_builder.hpp"

using namespace std;
using namespace emdw;

class LoopyAssocTest : public testing::Test {
	protected:
		virtual void SetUp() {
			inplaceNormalizer_ = uniqptr<FactorOperator> (new DiscreteTable_InplaceNormalize<T>);
			normalizer_ = uniqptr<FactorOperator> (new DiscreteTable_MaxNormalize<T>);
			marginalizer_ = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);

			// Probability table for a0
			a0Dom_ = uniqptr<DASS>(new DASS{0, 1});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{1}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a0}, {a0Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, marginalizer_, inplaceNormalizer_, normalizer_) ) );

			// Probability table for a1
			a1Dom_ = uniqptr<DASS>(new DASS{0, 1, 2});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{1}] = sparseProbs_[DASS{2}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a1}, {a1Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, marginalizer_, inplaceNormalizer_, normalizer_) ) );

			// Probability table for a2
			a2Dom_ = uniqptr<DASS>(new DASS{0, 2});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{2}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a2}, {a2Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, marginalizer_, inplaceNormalizer_, normalizer_) ) );

		}

		virtual void TearDown() {
			hypotheses_.clear();
			sparseProbs_.clear();
		}

	protected:
		typedef unsigned short T;
		typedef DiscreteTable<T> DT;
		typedef std::vector<T> DASS;

	protected:
		enum{a0, a1, a2};
		const double kFloor_ = 0.0;
		const double kMargin_ = 0.0;
		const double kDefProb_ = 0.0;

		rcptr<FactorOperator> inplaceNormalizer_;
		rcptr<FactorOperator> normalizer_;
		rcptr<FactorOperator> marginalizer_;

		rcptr<DASS> a0Dom_;
		rcptr<DASS> a1Dom_;
		rcptr<DASS> a2Dom_;

		vector<rcptr<Factor>> hypotheses_;
		map<DASS, FProb> sparseProbs_;
}; // LoopyAssocTest


// GraphBuilder tests
TEST_F (LoopyAssocTest, SmallTest) {
	// Association hypotheses
	std::map<RVIdType, rcptr<DASS>> assocHypotheses;
	assocHypotheses[1] = uniqptr<DASS>(new DASS{0, 1});
	assocHypotheses[2] = uniqptr<DASS>(new DASS{0, 2});
	assocHypotheses[3] = uniqptr<DASS>(new DASS{0, 1, 2});

	// Build the graphs
	rcptr<GraphBuilder> gb = uniqptr<GraphBuilder> (new GraphBuilder());
	std::vector<rcptr<Graph>> graph = gb->getGraphs( assocHypotheses );
	unsigned i = 1;

	// Marginal over var i
	//graph[0]->plotGraph();
	rcptr<Factor> prevMarginal = graph[0]->getMarginalBelief(i);

	// Message passing
	//for(unsigned i = 0; i < 25; i++) graph[0]->depthFirstMessagePassing();

	// Marginal over var i
	//graph[0]->plotGraph();
	rcptr<Factor> postMarginal = graph[0]->getMarginalBelief(i);

	EXPECT_EQ(0, 0);
} // GraphBuilderInit()

TEST_F (LoopyAssocTest, FullJointSanityTest) {
	rcptr<Factor> table;
	emdw::RVIds vars = {1, 2, 3, 4, 5};

	// Association hypotheses
	std::map<RVIdType, rcptr<DASS>> assocHypotheses;
	assocHypotheses[1] = uniqptr<DASS>(new DASS{0, 1});
	assocHypotheses[2] = uniqptr<DASS>(new DASS{0, 2});
	assocHypotheses[3] = uniqptr<DASS>(new DASS{0, 3});
	assocHypotheses[4] = uniqptr<DASS>(new DASS{0, 1, 3});
	assocHypotheses[5] = uniqptr<DASS>(new DASS{0, 1, 2, 3});

	// Create factors
	std::vector<rcptr<Factor>> factors(5);;
	for (unsigned i = 0; i < 5; i++) {
		std::map<DASS, FProb> sparseProbs;
		rcptr<DASS> aDom = assocHypotheses[vars[i]];

		sparseProbs[ DASS{ (*aDom)[0] } ] = 0.5;
		for (unsigned j = 1; j < aDom->size(); j++ ) sparseProbs[DASS{ (*aDom)[j]  }] = 1;
			
		factors[i] = uniqptr<Factor> (new DT(emdw::RVIds{vars[i]}, {aDom}, kDefProb_,
					sparseProbs, kMargin_, kFloor_, false,
					marginalizer_, inplaceNormalizer_, normalizer_ ) );

		sparseProbs.clear();
	}

	// Create full joint
	std::vector<rcptr<Factor>> pairwise(6);

	// Pairwise factor (1, 4)
	pairwise[0] = factors[0]->absorb(factors[3]);
	rcptr<DT> cancel = std::dynamic_pointer_cast<DT>(pairwise[0]);
	cancel->setEntry(emdw::RVIds{1, 4}, emdw::RVVals{T(1), T(1)}, 0);
	pairwise[0]->inplaceNormalize();

	// Pairwise factor (1, 5)
	pairwise[1] = factors[0]->absorb(factors[4]);
	cancel = std::dynamic_pointer_cast<DT>(pairwise[1]);
	cancel->setEntry(emdw::RVIds{1, 5}, emdw::RVVals{T(1), T(1)}, 0);
	pairwise[1]->inplaceNormalize();

	// Pairwise factor (2, 5)
	pairwise[2] = factors[1]->absorb(factors[4]);
	cancel = std::dynamic_pointer_cast<DT>(pairwise[2]);
	cancel->setEntry(emdw::RVIds{2, 5}, emdw::RVVals{T(2), T(2)}, 0);
	pairwise[2]->inplaceNormalize();

	// Pairwise factor (3, 4)
	pairwise[3] = factors[2]->absorb(factors[3]);
	cancel = std::dynamic_pointer_cast<DT>(pairwise[3]);
	cancel->setEntry(emdw::RVIds{3, 4}, emdw::RVVals{T(3), T(3)}, 0);
	pairwise[3]->inplaceNormalize();
	
	// Pairwise factor (3, 5)
	pairwise[4] = factors[2]->absorb(factors[4]);
	cancel = std::dynamic_pointer_cast<DT>(pairwise[4]);
	cancel->setEntry(emdw::RVIds{3, 5}, emdw::RVVals{T(3), T(3)}, 0);
	pairwise[4]->inplaceNormalize();

	// Pairwise factor (4, 5)
	pairwise[5] = factors[3]->absorb(factors[4]);
	cancel = std::dynamic_pointer_cast<DT>(pairwise[5]);
	cancel->setEntry(emdw::RVIds{4, 5}, emdw::RVVals{T(1), T(1)}, 0);
	cancel->setEntry(emdw::RVIds{4, 5}, emdw::RVVals{T(3), T(3)}, 0);
	pairwise[5]->inplaceNormalize();

	rcptr<ClusterGraph> clusterGraph = uniqptr<ClusterGraph>( new ClusterGraph( pairwise ) );
	//clusterGraph->exportToGraphViz("pairwise");

	std::map<Idx2, rcptr<Factor>> msgs; msgs.clear();
	MessageQueue msgQ; msgQ.clear();

	unsigned nMsgs = loopyBU_CG(*clusterGraph, msgs, msgQ, 0.0);
	for (unsigned i = 0; i < clusterGraph->factorPtrs.size(); i++ ) {
		std::cout << "The distribution over message " << i << " is "
			<< *queryLBU_CG(*clusterGraph, msgs, clusterGraph->factorPtrs[i]->getVars() ) << std::endl;
	}

	rcptr<Factor> joint = queryLBU_CG(*clusterGraph, msgs, clusterGraph->factorPtrs[4]->getVars() );
	//rcptr<Factor> marginal = (joint->marginalize(emdw::RVIds{3}))->normalize(); 
	//std::cout << *joint << std::endl;

	EXPECT_EQ(0, 0);
}

TEST_F (LoopyAssocTest, LargeTest) {
	// Association hypotheses
	std::map<RVIdType, rcptr<DASS>> assocHypotheses;
	assocHypotheses[1] = uniqptr<DASS>(new DASS{0, 1});
	assocHypotheses[2] = uniqptr<DASS>(new DASS{0, 2});
	assocHypotheses[3] = uniqptr<DASS>(new DASS{0, 3});
	assocHypotheses[4] = uniqptr<DASS>(new DASS{0, 1, 3});
	assocHypotheses[5] = uniqptr<DASS>(new DASS{0, 1, 2, 3});
	assocHypotheses[6] = uniqptr<DASS>(new DASS{0, 4, 5});
	assocHypotheses[7] = uniqptr<DASS>(new DASS{0, 4});

	// Build the graphs
	rcptr<GraphBuilder> gb = uniqptr<GraphBuilder> (new GraphBuilder());
	std::vector<rcptr<Graph>> graph = gb->getGraphs( assocHypotheses );
	unsigned i = 3;

	// Marginal belief over var i
	//graph[0]->plotGraph();
	rcptr<Factor> prevMarginal = graph[0]->getMarginalBelief(i);
	//std::cout << *prevMarginal << std::endl;

	// Message passing
	for (unsigned i = 0; i < 1; i++) graph[0]->depthFirstMessagePassing();

	//graph[0]->plotGraph();
	
	// Marginal belief over var i
	rcptr<Factor> postMarginal = graph[0]->getMarginalBelief(i);
	//std::cout << *postMarginal << std::endl;

	EXPECT_EQ(0, 0);
} // GraphBuilderInit()

// Explicit handmade example, just to get a rough idea of the process.
TEST_F (LoopyAssocTest, HandHolding) {	
	vector<rcptr<Factor>> clusters;
	clusters.push_back(hypotheses_[0]->absorb(hypotheses_[1]));
	clusters.push_back(hypotheses_[0]->absorb(hypotheses_[2]));
	clusters.push_back(hypotheses_[1]->absorb(hypotheses_[2]));

	enum{c0, c1, c2};
	std::map<RVIds, rcptr<Factor>> old_messages;
	old_messages.clear();

	old_messages[RVIds{c0, c1}] = rcptr<Factor> (clusters[0]->vacuousCopy(RVIds{a0}));
	old_messages[RVIds{c1, c2}] = rcptr<Factor> (clusters[1]->vacuousCopy(RVIds{a2}));
	old_messages[RVIds{c2, c0}] = rcptr<Factor> (clusters[2]->vacuousCopy(RVIds{a1}));


	//for (unsigned i = 0; i < 3; i++) cout << *clusters[i] << endl;

	// Cluster 0 generates message 1
	rcptr<DT> change = std::dynamic_pointer_cast<DT>(clusters[0]);
	change->setEntry(RVIds{a0, a1}, RVVals{T(0), T(0)}, 0);
	change->setEntry(RVIds{a0, a1}, RVVals{T(1), T(1)}, 0);
	rcptr<Factor> msg = clusters[0]->marginalize(RVIds{a0});
	
	//cout << "a0 marginal" << endl;
	//cout << *msg << endl;
	
	// Cluster 1 absorbs message 1
	clusters[1]->inplaceAbsorb(msg);
	change = std::dynamic_pointer_cast<DT>(clusters[1]);
	change->setEntry(RVIds{a0, a2}, RVVals{T(0), T(0)}, 0);
	msg = clusters[1]->marginalize(RVIds{a2});

	//cout << "a2 marginal" << endl;
	//cout << *msg << endl;
	
	// Cluster 2 absorbs message 2
	clusters[2]->inplaceAbsorb(msg);
	change = std::dynamic_pointer_cast<DT>(clusters[2]);
	change->setEntry(RVIds{a1, a2}, RVVals{T(0), T(0)}, 0);
	change->setEntry(RVIds{a1, a2}, RVVals{T(2), T(2)}, 0);
	msg = clusters[2]->marginalize(RVIds{a1});

	//cout << "a1 marginal" << endl;
	//cout << *msg << endl;

	//cout << "a2 marginal" << endl;
	//cout << *clusters[2]->marginalize(RVIds{a2}) << endl;

	// Cluster 0 absorbs message 3
	
	clusters[0]->inplaceAbsorb(msg);
	change = std::dynamic_pointer_cast<DT>(clusters[0]);
	change->setEntry(RVIds{a0, a1}, RVVals{T(0), T(0)}, 0);
	change->setEntry(RVIds{a0, a1}, RVVals{T(1), T(1)}, 0);
	msg = clusters[0]->marginalize(RVIds{a0});

	//cout << "a0 marginal" << endl;
	//cout << *msg << endl;

	EXPECT_EQ(0, 0);
} // HandHolding()
