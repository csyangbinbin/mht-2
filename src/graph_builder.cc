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
			ConstructClusters();
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
		
		distribution_[a_[i]] = uniqptr<DT> (new DT(emdw::RVIds{a_[i]}, {a_dom}, def_prob_,
						sparse_probs, margin_, floor_, false,
						marg_ptr_, inorm_ptr_, norm_ptr_ ) );
		
		sparse_probs.clear();
	}
} // ConstructDistributions()

void GraphBuilder::ConstructClusters() {
	DASS intersection;
	emdw::RVIds neighbours;

	std::vector<emdw::RVIds> neighbourhoods(a_.size());
	std::vector<emdw::RVIds> pairs;
	
	rcptr<DASS> ai_dom;
	rcptr<DASS> aj_dom;

	// Step 1: Group all association hypotheses which share common origins.
	for (unsigned i = 0; i < a_.size(); i++) {
		
		ai_dom = assoc_hypotheses_[a_[i]]; 
		
		for (unsigned j = 0; j < a_.size(); j++) {
			aj_dom = assoc_hypotheses_[a_[j]];

			std::set_intersection((*ai_dom).begin(), (*ai_dom).end(),
					(*aj_dom).begin(), (*aj_dom).end(),
					std::back_inserter(intersection));

			if (intersection.size() > 1) neighbourhoods[i].push_back(a_[j]);
			intersection.clear(); 
		}
	}

	// Step 2: Generate all pairwise connections implied by shared associations
	for (unsigned i = 0; i < neighbourhoods.size(); i++) {
		neighbours = neighbourhoods[i];
		for (unsigned j = 0; j < neighbours.size(); j++) {
			for (unsigned k = j+1; k < neighbours.size(); k++) {
				pairs.push_back(emdw::RVIds{neighbours[j], neighbours[k]});
			}
			if (neighbours.size() == 1) pairs.push_back(emdw::RVIds{neighbours[j], neighbours[j]});
		}
	}

	// Step 3: Select only the unique pairings
	std::sort(pairs.begin(), pairs.end());
	pairs.erase(unique(pairs.begin(), pairs.end()), pairs.end());

	//for (unsigned i = 0; i < pairs.size(); i++) std::cout << pairs[i] << std::endl;
	
	// Step 4: Create all pairwise (or single) clusters
	for (unsigned i = 0; i < pairs.size(); i++) {
		cluster_.push_back( distribution_[pairs[i][0]]->absorb(distribution_[pairs[i][1]]) );
	}

	for (unsigned i = 0; i < cluster_.size(); i++) std::cout << *cluster_[i] << std::endl;

} // ConstructClusters()
