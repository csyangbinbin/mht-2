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
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "transforms.hpp"
#include "utils.hpp"

using namespace std;
using namespace emdw;

TEST (LoopyAssocTest, Loopy) {
	// Variable type declaration
	typedef unsigned short T;
	typedef DiscreteTable<T> DT;
	typedef vector<T> DASS;
	
	double floor = 0.0;
	double margin = 2.0*floor;
	double def_prob = 0.0;

	// Factor Operator declaration
	rcptr<FactorOperator> marg_ptr = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);
	rcptr<FactorOperator> inorm_ptr = uniqptr<FactorOperator> (new DiscreteTable_InplaceNormalize<T>);
	rcptr<FactorOperator> norm_ptr = uniqptr<FactorOperator> (new DiscreteTable_MaxNormalize<T>);

	// Random variable name declaration
	enum{a0, a1, a2};
	rcptr<DASS> a0_dom (new DASS{0, 1});
	vector<rcptr<Factor>> hypotheses;

	// Probability table for a0
	map<DASS, FProb> sparseProbs; sparseProbs.clear();
	sparseProbs[DASS{0}] = sparseProbs[DASS{1}] = 1;

	// Discrete probability table
	hypotheses.push_back(uniqptr<DT> (new DT(RVIds{a0}, {a0_dom}, def_prob, 
					sparseProbs, margin, floor, false, marg_ptr, inorm_ptr, norm_ptr) ) );
	
	
	
	// Probability table for a1
	rcptr<DASS> a1_dom (new DASS{0, 1, 2});
	sparseProbs.clear();
	sparseProbs[DASS{0}] = sparseProbs[DASS{1}] = sparseProbs[DASS{2}] = 1;
	hypotheses.push_back(uniqptr<DT> (new DT(RVIds{a1}, {a1_dom}, def_prob, 
					sparseProbs, margin, floor, false, marg_ptr, inorm_ptr, norm_ptr) ) );


	// Probability table for a2
	rcptr<DASS> a2_dom (new DASS{0, 2});
	sparseProbs.clear(); 
	sparseProbs[DASS{0}] = sparseProbs[DASS{2}] = 1;
	hypotheses.push_back(uniqptr<DT> (new DT(RVIds{a2}, {a2_dom}, def_prob, 
					sparseProbs, margin, floor, false, marg_ptr, inorm_ptr, norm_ptr) ) );

	
	vector<rcptr<Factor>> clusters;
	clusters.push_back(hypotheses[0]->absorb(hypotheses[1]));
	clusters.push_back(hypotheses[0]->absorb(hypotheses[2]));
	clusters.push_back(hypotheses[1]->absorb(hypotheses[2]));

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
	
	cout << "a0 marginal" << endl;
	cout << *msg << endl;
	
	// Cluster 1 absorbs message 1
	clusters[1]->inplaceAbsorb(msg);
	change = std::dynamic_pointer_cast<DT>(clusters[1]);
	change->setEntry(RVIds{a0, a2}, RVVals{T(0), T(0)}, 0);
	msg = clusters[1]->marginalize(RVIds{a2});

	cout << "a2 marginal" << endl;
	cout << *msg << endl;
	
	// Cluster 2 absorbs message 2
	clusters[2]->inplaceAbsorb(msg);
	change = std::dynamic_pointer_cast<DT>(clusters[2]);
	change->setEntry(RVIds{a1, a2}, RVVals{T(0), T(0)}, 0);
	change->setEntry(RVIds{a1, a2}, RVVals{T(2), T(2)}, 0);
	msg = clusters[2]->marginalize(RVIds{a1});

	cout << "a1 marginal" << endl;
	cout << *msg << endl;

	cout << "a2 marginal" << endl;
	cout << *clusters[2]->marginalize(RVIds{a2}) << endl;

	// Cluster 0 absorbs message 3
	
	clusters[0]->inplaceAbsorb(msg);
	change = std::dynamic_pointer_cast<DT>(clusters[0]);
	change->setEntry(RVIds{a0, a1}, RVVals{T(0), T(0)}, 0);
	change->setEntry(RVIds{a0, a1}, RVVals{T(1), T(1)}, 0);
	msg = clusters[0]->marginalize(RVIds{a0});

	cout << "a0 marginal" << endl;
	cout << *msg << endl;

	EXPECT_EQ(0, 0);
}
