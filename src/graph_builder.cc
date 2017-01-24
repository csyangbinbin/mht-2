/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Implements the graph builder class. For now it creates a singly connected
 * (in terms of sepsets) cluster graph. The nodes are pairwise unnormalised
 * measures.
 *************************************************************************/

#include <vector>
#include <iostream>
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "matops.hpp"
#include "vecset.hpp"
#include "graph_builder.hpp"

GraphBuilder::GraphBuilder(const std::map<emdw::RVIdType, rcptr<DASS>>& assoc_hypotheses) :
		assoc_hypotheses_(assoc_hypotheses) {
			rv_ids_ = this->getRVIds();
			constructDistributions();
} // GraphBuilder()

GraphBuilder::GraphBuilder() {} // GraphBuilder()

GraphBuilder::~GraphBuilder() {} // ~GraphBuilder()

emdw::RVIds GraphBuilder::getRVIds() const {
	emdw::RVIds a;

	for (auto const& e : this->assoc_hypotheses_) {
		a.push_back(e.first);
		std::sort((e.second)->begin(), (e.second)->end());
	}
	
	return a;
} // getRVIds()

void GraphBuilder::constructDistributions() {
	emdw::RVIds a = this->rv_ids_;
	std::vector <rcptr<Factor>> distributions;
	std::map<DASS, FProb> sparse_probs;

	rcptr<FactorOperator> marg_ptr = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);
	rcptr<FactorOperator> inorm_ptr = uniqptr<FactorOperator> (new DiscreteTable_InplaceNormalize<T>);
	rcptr<FactorOperator> norm_ptr = uniqptr<FactorOperator> (new DiscreteTable_Normalize<T>);

	double floor = 0;
	double margin = 0;
	double def_prob = 0;

	for (unsigned i = 0; i < a.size(); i++) {
		rcptr<DASS> a_dom = assoc_hypotheses_[a[i]];
		
		for (unsigned j = 0; j < a_dom->size(); j++) sparse_probs[DASS{(*a_dom)[j]}] = 1;
		
		distributions.push_back(uniqptr<DT> (new DT(emdw::RVIds{a[i]}, {a_dom}, def_prob,
						sparse_probs, margin, floor, false,
						marg_ptr, inorm_ptr, norm_ptr ) ) );
		
		sparse_probs.clear();
	}

} // constructDistributions()

void GraphBuilder::constructClusters() {
	emdw::RVIds a = this->rv_ids_;

	for (unsigned i = 0; i < a.size(); i++) {
		for (unsigned j = i+1; j < a.size(); j++) {
			std::cout << i << " - " << j << std::endl;
		}
	}

} // constructClusters()
