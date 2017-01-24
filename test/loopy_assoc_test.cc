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
#include "graph_builder.hpp"

using namespace std;
using namespace emdw;

class LoopyAssocTest : public testing::Test {
	protected:
		virtual void SetUp() {
			marg_ptr_ = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);
			inorm_ptr_ = uniqptr<FactorOperator> (new DiscreteTable_InplaceNormalize<T>);
			norm_ptr_ = uniqptr<FactorOperator> (new DiscreteTable_MaxNormalize<T>);

			// Probability table for a0
			a0_dom_ = uniqptr<DASS>(new DASS{0, 1});
			sparse_probs_.clear();
			sparse_probs_[DASS{0}] = sparse_probs_[DASS{1}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a0}, {a0_dom_}, kDefProb_, 
					sparse_probs_, kMargin_, kFloor_, false, marg_ptr_, inorm_ptr_, norm_ptr_) ) );

			// Probability table for a1
			a1_dom_ = uniqptr<DASS>(new DASS{0, 1, 2});
			sparse_probs_.clear();
			sparse_probs_[DASS{0}] = sparse_probs_[DASS{1}] = sparse_probs_[DASS{2}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a1}, {a1_dom_}, kDefProb_, 
					sparse_probs_, kMargin_, kFloor_, false, marg_ptr_, inorm_ptr_, norm_ptr_) I) );

			// Probability table for a2
			a2_dom_ = uniqptr<DASS>(new DASS{0, 2});
			sparse_probs_.clear();
			sparse_probs_[DASS{0}] = sparse_probs_[DASS{2}] = 1;

			hypotheses_.push_back(uniqptr<DT> (new DT(RVIds{a2}, {a2_dom_}, kDefProb_, 
					sparse_probs_, kMargin_, kFloor_, false, marg_ptr_, inorm_ptr_, norm_ptr_) ) );

		}

		virtual void TearDown() {
			hypotheses_.clear();
			sparse_probs_.clear();
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

		rcptr<FactorOperator> marg_ptr_;
		rcptr<FactorOperator> inorm_ptr_;
		rcptr<FactorOperator> norm_ptr_;

		rcptr<DASS> a0_dom_;
		rcptr<DASS> a1_dom_;
		rcptr<DASS> a2_dom_;

		vector<rcptr<Factor>> hypotheses_;
		map<DASS, FProb> sparse_probs_;
}; // LoopyAssocTest


TEST_F (LoopyAssocTest, GraphBuilder) {

	std::map<RVIdType, rcptr<DASS>> assoc_hypotheses;
	assoc_hypotheses[a0] = a0_dom_;
	assoc_hypotheses[a1] = a1_dom_;
	assoc_hypotheses[a2] = a2_dom_;

	rcptr<GraphBuilder> gb = uniqptr<GraphBuilder> (new GraphBuilder(assoc_hypotheses));
	
	EXPECT_EQ(0, 0);
}


 // Explicit handmade example, just to get a rough idea of the process.
TEST_F (LoopyAssocTest, Loopy) {	
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
}
