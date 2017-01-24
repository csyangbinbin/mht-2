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

GraphBuilder::GraphBuilder(const std::map<emdw::RVIdType, rcptr<DASS>>& assoc_hypotheses,
		const rcptr<FactorOperator>& marg_ptr, 
		const rcptr<FactorOperator>& inorm_ptr,
		const rcptr<FactorOperator>& norm_ptr,
		const double floor,
		const double margin,
		const double def_prob) :
		assoc_hypotheses_(assoc_hypotheses),
		marg_ptr_(marg_ptr),
		inorm_ptr_(inorm_ptr),
		norm_ptr_(norm_ptr),
		floor_(floor),
		margin_(margin),
		def_prob_(def_prob) {
			AcquireRVIds();
			ConstructDistributions();
} // GraphBuilder()

GraphBuilder::GraphBuilder() {} // GraphBuilder()

GraphBuilder::~GraphBuilder() {} // ~GraphBuilder()

void GraphBuilder::AcquireRVIds() {
	for (auto const& e : assoc_hypotheses_) {
		a_.push_back(e.first);
		std::sort((e.second)->begin(), (e.second)->end());
	}
} // AcquireRVIds()

void GraphBuilder::ConstructDistributions() {
	std::map<DASS, FProb> sparse_probs;

	for (unsigned i = 0; i < a_.size(); i++) {
		rcptr<DASS> a_dom = assoc_hypotheses_[a_[i]];
		
		for (unsigned j = 0; j < a_dom->size(); j++) sparse_probs[DASS{(*a_dom)[j]}] = 1;
		
		distribution_.push_back(uniqptr<DT> (new DT(emdw::RVIds{a_[i]}, {a_dom}, def_prob_,
						sparse_probs, margin_, floor_, false,
						marg_ptr_, inorm_ptr_, norm_ptr_ ) ) );
		
		sparse_probs.clear();
	}
} // ConstructDistributions()

void GraphBuilder::ConstructClusters() {
	for (unsigned i = 0; i < a_.size(); i++) {
		for (unsigned j = i+1; j < a_.size(); j++) {
			std::cout << i << " - " << j << std::endl;
		}
	}
} // ConstructClusters()
