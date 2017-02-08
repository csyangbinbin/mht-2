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


rcptr<FactorOperator> defaultInplaceNormalizerGB = uniqptr<FactorOperator>(new DiscreteTable_InplaceNormalize<unsigned short>);
rcptr<FactorOperator> defaultNormalizerGB = uniqptr<FactorOperator>(new DiscreteTable_MaxNormalize<unsigned short>);
rcptr<FactorOperator> defaultMarginalizerGB = uniqptr<FactorOperator>(new DiscreteTable_Marginalize<unsigned short>);

GraphBuilder::GraphBuilder(
		const std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses,
		const double floor,
		const double margin,
		const double defProb,
		const rcptr<FactorOperator>& inplaceNormalizer, 
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& marginalizer) 
			: assocHypotheses_(assocHypotheses),
			  floor_(floor),
			  margin_(margin),
			  defProb_(defProb) {
			
	// Default initialisation
	if(!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerGB;
	if(!normalizer) normalizer_ = defaultNormalizerGB;
	if(!marginalizer) marginalizer_ = defaultMarginalizerGB;

	extractRVIds();
	constructDistributions();
	constructClusters();

} // GraphBuilder()

GraphBuilder::GraphBuilder() {} // GraphBuilder()

GraphBuilder::~GraphBuilder() {} // ~GraphBuilder()

void GraphBuilder::extractRVIds() {
	for (auto const& e : assocHypotheses_) {
		a_.push_back(e.first);
		std::sort((e.second)->begin(), (e.second)->end());
	}
} // AcquireRVIds()

void GraphBuilder::constructDistributions() {
	std::map<DASS, FProb> sparseProbs;
	rcptr<DASS> aDom;

	for (unsigned i = 0; i < a_.size(); i++) {
		aDom = assocHypotheses_[a_[i]];
		
		for (unsigned j = 0; j < aDom->size(); j++) sparseProbs[DASS{(*aDom)[j]}] = 1;
		
		dist_[a_[i]] = uniqptr<DT> (new DT(emdw::RVIds{a_[i]}, {aDom}, defProb_,
						sparseProbs, margin_, floor_, false,
						marginalizer_, inplaceNormalizer_, normalizer_ ) );
		
		sparseProbs.clear();
	}
} // ConstructDistributions()

void GraphBuilder::constructClusters() {
	DASS domIntersection;
	emdw::RVIds nodeVars, maximalClique, varIntersection;

	std::vector<rcptr<Node>> vertex;
	std::vector<emdw::RVIds> clique(a_.size()), districts;
	rcptr<DASS> aiDom, ajDom;

	bool intersects = false;
	int prev = -1;

	// Step 1: Group all association hypotheses which share common origins.
	for (unsigned i = 0; i < a_.size(); i++) {
		
		aiDom = assocHypotheses_[a_[i]]; 
		
		for (unsigned j = 0; j < a_.size(); j++) {
			ajDom = assocHypotheses_[a_[j]];

			std::set_intersection(aiDom->begin(), aiDom->end(),
					ajDom->begin(), ajDom->end(),
					std::back_inserter(domIntersection));

			if (domIntersection.size() > 1) clique[i].push_back(a_[j]);
			domIntersection.clear(); 
		}
	}

	// Step 2: Determine all disjoint association variable networks
	for (unsigned i = 0; i < clique.size(); i++) {
		for (unsigned j = 0; j < clique.size(); j++) {

			std::set_intersection( clique[i].begin(), clique[j].end(),
					clique[j].begin(), clique[j].end(),
					std::back_inserter(varIntersection));

			if (varIntersection.size() > 0) {
				std::set_union( clique[i].begin(), clique[i].end(),
						clique[j].begin(), clique[j].end(),
						std::back_inserter(maximalClique));

				clique[i] = std::move(maximalClique);
				maximalClique.clear();
			}
			varIntersection.clear();
		}
	}

	// Step 3: Select only unique cliques
	std::sort(clique.begin(), clique.end());
	clique.erase(unique( clique.begin(), clique.end() ), clique.end());

	// Step 4: Generate all pairwise connections in the clique
	std::vector<std::vector<emdw::RVIds>> pairs(clique.size());
	unsigned numNodes;
	
	for (unsigned i = 0; i < clique.size(); i++) {	
		
		nodeVars = clique[i]; numNodes = nodeVars.size();

		for (unsigned j = 0; j < numNodes; j++) {
			for (unsigned k = j+1; k < numNodes; k++) pairs[i].push_back(emdw::RVIds{nodeVars[j], nodeVars[k]});
			if (numNodes == 1) pairs[i].push_back(emdw::RVIds{nodeVars[j], nodeVars[j]});
		}
	}
	
	// Step 5: Select only the unique pairings
	for (auto& p : pairs) {
		std::sort(p.begin(), p.end());
		p.erase( unique(p.begin(), p.end()), p.end() );
	}

	// Step 6: Create the pairwise clusters
	graphs_.clear(); graphs_.resize(pairs.size());

	for (unsigned i = 0; i < pairs.size(); i++) {
		graphs_[i] = uniqptr<Graph>(new Graph());

		// Step 6.1 : Create the nodes for the graph
		for (auto& j : pairs[i]) vertex.push_back(uniqptr<Node>(new Node ( ( dist_[j[0]] )->absorb( dist_[j[1]] ) ) ) );

		// Step 6.2 : Link up the nodes
		for (auto& j : clique[i]) {
			for (unsigned k = 0; k < vertex.size(); k++) {
   				nodeVars = vertex[k]->getVars();
				
				// If the cluster contains variable j
				for (auto& v : nodeVars)  {
					if (v == j) { intersects = true; }
				}

				// Link up the nodes in a chain
				if (intersects && prev != -1) graphs_[i]->addEdge(vertex[prev], vertex[k]); 
				if (intersects) prev = (int) k;

				// Clear the values
				intersects = false;
				nodeVars.clear();
			}
			prev = -1;
		}
		vertex.clear();	
	}

	graphs_[1]->depthFirstSearch();

} // ConstructClusters()
