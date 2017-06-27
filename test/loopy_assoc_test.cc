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
			normalizer_ = uniqptr<FactorOperator> (new DiscreteTable_Normalize<T>);
			marginalizer_ = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);
		}

		virtual void TearDown() {}

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
}; // LoopyAssocTest


// GraphBuilder tests
TEST_F (LoopyAssocTest, SmallTest) {
	emdw::RVIds vars = {1, 2, 3};

	// Association hypotheses
	std::map<RVIdType, rcptr<DASS>> assocHypotheses;
	assocHypotheses[1] = uniqptr<DASS>(new DASS{0, 1});
	assocHypotheses[2] = uniqptr<DASS>(new DASS{0, 2});
	assocHypotheses[3] = uniqptr<DASS>(new DASS{0, 1, 2});

	// Build the graphs
	rcptr<GraphBuilder> gb = uniqptr<GraphBuilder> (new GraphBuilder());
	std::map<emdw::RVIdType, rcptr<Factor>> marginals = gb->getMarginals( assocHypotheses );

	assocHypotheses[1] = uniqptr<DASS>(new DASS{0, 1});
	assocHypotheses[2] = uniqptr<DASS>(new DASS{0, 2});
	/*
	for (emdw::RVIdType i : vars) {
		 std::cout << "Belief held over var " << i << "\n" << *(marginals[i]) << 
			 "=====================================" << std::endl;
	}
	*/

	EXPECT_EQ(0, 0);
} // SmallTest()

TEST_F (LoopyAssocTest, LargeTest) {
	emdw::RVIds vars = {1, 2, 3, 4, 5};

	// Association hypotheses
	std::map<RVIdType, rcptr<DASS>> assocHypotheses;
	assocHypotheses[1] = uniqptr<DASS>(new DASS{0, 1});
	assocHypotheses[2] = uniqptr<DASS>(new DASS{0, 3});
	assocHypotheses[3] = uniqptr<DASS>(new DASS{0, 2});
	assocHypotheses[4] = uniqptr<DASS>(new DASS{0, 1, 2});
	assocHypotheses[5] = uniqptr<DASS>(new DASS{0, 1, 2, 3});
	//assocHypotheses[6] = uniqptr<DASS>(new DASS{0, 4, 5});
	//assocHypotheses[7] = uniqptr<DASS>(new DASS{0, 4});

	// Build the graphs
	rcptr<GraphBuilder> gb = uniqptr<GraphBuilder> (new GraphBuilder());
	std::map<emdw::RVIdType, rcptr<Factor>> marginals = gb->getMarginals( assocHypotheses );

	/*
	for (emdw::RVIdType i : vars) {
		 std::cout << "Belief held over var " << i << "\n" << *(marginals[i]) << 
			 "=====================================" << std::endl;
	}
	*/

	EXPECT_EQ(0, 0);
} // LargeTest()

TEST_F (LoopyAssocTest, FullJointSanityTest) {
	rcptr<Factor> table;
	emdw::RVIds vars = {1, 2, 3, 4, 5, 6, 7};

	// Association hypotheses
	std::map<RVIdType, rcptr<DASS>> assocHypotheses;
	assocHypotheses[1] = uniqptr<DASS>(new DASS{0, 1});
	assocHypotheses[2] = uniqptr<DASS>(new DASS{0, 2});
	assocHypotheses[3] = uniqptr<DASS>(new DASS{0, 3});
	assocHypotheses[4] = uniqptr<DASS>(new DASS{0, 1, 3});
	assocHypotheses[5] = uniqptr<DASS>(new DASS{0, 1, 2, 3});
	assocHypotheses[6] = uniqptr<DASS>(new DASS{0, 4, 5});
	assocHypotheses[7] = uniqptr<DASS>(new DASS{0, 4});

	// Create factors
	std::vector<rcptr<Factor>> factors(7);
	for (unsigned i = 0; i < 7; i++) {
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
	std::vector<rcptr<Factor>> pairwise(7);

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

	// Pairwise factor (4, 5)
	pairwise[6] = factors[5]->absorb(factors[6]);
	cancel = std::dynamic_pointer_cast<DT>(pairwise[6]);
	cancel->setEntry(emdw::RVIds{6, 7}, emdw::RVVals{T(4), T(4)}, 0);
	pairwise[6]->inplaceNormalize();

	rcptr<ClusterGraph> clusterGraph = uniqptr<ClusterGraph>( new ClusterGraph( pairwise ) );
	//clusterGraph->exportToGraphViz("pairwise");

	std::map<Idx2, rcptr<Factor>> msgs; msgs.clear();
	MessageQueue msgQ; msgQ.clear();

	unsigned nMsgs = loopyBU_CG(*clusterGraph, msgs, msgQ, 0.0);
	
	/*
	for (emdw::RVIdType i : vars) {
		 std::cout << "Belief held over var " << i << "\n" <<
		 *(queryLBU_CG(*clusterGraph, msgs, emdw::RVIds{ i  } )->normalize()) <<
			 "=====================================" << std::endl;
	}
	*/

	EXPECT_EQ(0, 0);
}
