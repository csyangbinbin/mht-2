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

//TODO: Add all the default DiscreteTable factor operators
rcptr<FactorOperator> defaultInplaceNormalizerGB = uniqptr<FactorOperator>(new DiscreteTable_InplaceNormalize<unsigned short>);
rcptr<FactorOperator> defaultNormalizerGB = uniqptr<FactorOperator>(new DiscreteTable_MaxNormalize<unsigned short>);
rcptr<FactorOperator> defaultMarginalizerGB = uniqptr<FactorOperator>(new DiscreteTable_Marginalize<unsigned short>);

GraphBuilder::GraphBuilder(
		const double floor,
		const double margin,
		const double defProb,
		const rcptr<FactorOperator>& inplaceNormalizer, 
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& marginalizer) 
			: floor_(floor),
			  margin_(margin),
			  defProb_(defProb),
			  inplaceNormalizer_(inplaceNormalizer),
			  normalizer_(normalizer),
			  marginalizer_(marginalizer)
	{			
	// Default initialisation
	if(!inplaceNormalizer_) inplaceNormalizer_ = defaultInplaceNormalizerGB;
	if(!normalizer_) normalizer_ = defaultNormalizerGB;
	if(!marginalizer_) marginalizer_ = defaultMarginalizerGB;
} // GraphBuilder()

GraphBuilder::~GraphBuilder() {} // Default Destructor

std::vector<rcptr<Graph>> GraphBuilder::getGraphs(std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const {
	std::vector<rcptr<Graph>> graphs;

	// Extract RVIds
	emdw::RVIds vars = extractRVIds(assocHypotheses);

	// Construct the distributions
	std::map<emdw::RVIdType, rcptr<Factor>> dist = constructDistributions(vars, assocHypotheses);
	
	// Construct the graphs
	graphs = constructClusters();

	return graphs;
} // getGraphs()


emdw::RVIds GraphBuilder::extractRVIds(const std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const {
	emdw::RVIds vars; vars.clear();

	for (auto const& e : assocHypotheses) {
		vars.push_back(e.first);
		std::sort((e.second)->begin(), (e.second)->end());
	}

	return vars;
} // extractRVIds()

std::map<emdw::RVIdType, rcptr<Factor>> GraphBuilder::constructDistributions(
		const emdw::RVIds& vars, 
		std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const {
	std::map<emdw::RVIdType, rcptr<Factor>> dist;

	for (unsigned i = 0; i < vars.size(); i++) {
		std::map<DASS, FProb> sparseProbs;
		rcptr<DASS> aDom = assocHypotheses[vars[i]];
		
		for (unsigned j = 0; j < aDom->size(); j++) sparseProbs[DASS{(*aDom)[j]}] = 1;
		
		dist[vars[i]] = uniqptr<DT> (new DT(emdw::RVIds{vars[i]}, {aDom}, defProb_,
						sparseProbs, margin_, floor_, false,
						marginalizer_, inplaceNormalizer_, normalizer_ ) );
		
		sparseProbs.clear();
	}

	return dist;
} // constructDistributions()

std::vector<rcptr<Graph>> GraphBuilder::constructClusters() const {
	std::vector<rcptr<Graph>> graphs(1);
	return graphs;
} // ConstructClusters()

