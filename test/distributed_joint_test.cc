/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for distributed joint test.
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
#include "combinations.hpp"
#include "graph_builder.hpp"

using namespace std;
using namespace emdw;

class DistributedJointTest : public testing::Test {
	protected:
		virtual void SetUp() {
			inplaceNormalizer_ = uniqptr<FactorOperator> (new DiscreteTable_InplaceNormalize<T>);
			normalizer_ = uniqptr<FactorOperator> (new DiscreteTable_Normalize<T>);
			marginalizer_ = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);
		}

		virtual void TearDown() {}

	protected:
		typedef unsigned T;
		typedef DiscreteTable<T> DT;
		typedef std::vector<T> DASS;

	protected:
		// Maximum number of targets and measurements
		const unsigned kNumTargets_ = 5;
		const unsigned kNumMeasurements_ = 4;

		// Parameters for the DiscreteTable RV
		const double kFloor_ = 0.0;
		const double kMargin_ = 0.0;
		const double kDefProb_ = 0.0;

		// Factor operators for DiscreteTable
		rcptr<FactorOperator> inplaceNormalizer_;
		rcptr<FactorOperator> normalizer_;
		rcptr<FactorOperator> marginalizer_;
}; // LoopyAssocTest


// GraphBuilder tests
TEST_F (DistributedJointTest, RunTest) {

	// Declare the global domains for an infinite gate
	DASS globalDomain(kNumMeasurements_);
	for (unsigned i = 0; i < kNumMeasurements_; i++) globalDomain[i] = i;

	// Declare variable identities and their respective domains
	emdw::RVIds vars(kNumTargets_);
	std::vector< rcptr<DASS> > domains(kNumTargets_);
	for (unsigned i = 0; i < kNumTargets_; i++) {
		vars[i] = i;
		domains[i] = uniqptr<DASS>( new DASS(globalDomain) );
	}

	std::vector< std::vector<unsigned> > combinations = getAllCombinations(kNumTargets_, kNumMeasurements_);
	unsigned N = combinations.size();
	
	// for (unsigned i = 0; i < N; i++) std::cout << combinations[i] << std::endl;

	std::map<DASS, FProb> sparseProbs; sparseProbs.clear();
	for (unsigned i = 0; i < N; i++) {
		std::vector<unsigned> count(kNumMeasurements_, 0);
		bool repeated = false;

		for (unsigned j = 0; j < kNumTargets_; j++) {
			unsigned k = combinations[i][j];
			repeated = (k != 0 && ++count[k] > 1);
			if (repeated) break;
		}
		
		if (repeated) continue;

		sparseProbs[combinations[i]] = 1.0;
	} 

	rcptr<Factor> exactJoint= uniqptr<DT> (new DT( vars, domains,
				kDefProb_, sparseProbs,
				kFloor_, kMargin_,
				true, marginalizer_, 
				inplaceNormalizer_, normalizer_ ) );

	exactJoint->inplaceNormalize();
	rcptr<Factor> exactMarginal = exactJoint->marginalize({0});

	std::cout << "The exact belief for a0: " << std::endl;
	std::cout << *exactMarginal << std::endl;

	std::vector< std::vector<unsigned> > pairwiseScopes = getUniqueCombinations(2, kNumTargets_);
	unsigned M = pairwiseScopes.size();
	
	std::vector< std::vector<unsigned> > pairwiseDomains = getUniquePermutations(2, kNumMeasurements_);
	pairwiseDomains.push_back( {0, 0} );
	unsigned P = pairwiseDomains.size();

	std::vector<rcptr<Factor>> factors(M);

	for (unsigned i = 0; i < M; i++) {
		RVIds scope(pairwiseScopes[i]);

		std::vector< rcptr<DASS> > localDomain(2);
		localDomain[0] = localDomain[1] = uniqptr<DASS>( new DASS(globalDomain) ); 

		std::map<DASS, FProb> localProbs; localProbs.clear();
		for (unsigned j = 0; j < P; j++) localProbs[ pairwiseDomains[j] ] = 1.0;

		factors[i] = uniqptr<DT> (new DT( scope, localDomain,
				kDefProb_, localProbs,
				kFloor_, kMargin_,
				true, marginalizer_, 
				inplaceNormalizer_, normalizer_ ) );
	}

	for (unsigned i = 0; i < M; i++) std::cout << pairwiseScopes[i] << std::endl;
	
	rcptr<ClusterGraph> clusterGraph = uniqptr<ClusterGraph>(new ClusterGraph(factors));
	clusterGraph->exportToGraphViz("distributed");
	map<Idx2, rcptr<Factor>> msgs;
	MessageQueue msgQ;

	msgs.clear(); msgQ.clear();
	unsigned nMsgs = loopyBU_CG(*clusterGraph, msgs, msgQ, 0.0);
	cout << "Sent " << nMsgs << " before convergence." << std::endl;

	std::cout << "The final belief for a0: " << std::endl;
	std::cout << *queryLBU_CG(*clusterGraph, msgs, {0}) << std::endl;


} 
