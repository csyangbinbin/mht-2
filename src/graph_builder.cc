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
	graphs = constructClusters(vars, assocHypotheses, dist);

	return graphs;
} // getGraphs()


emdw::RVIds GraphBuilder::extractRVIds(const std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const {
	emdw::RVIds vars; vars.clear();

	for (auto const& e : assocHypotheses) {
		vars.push_back(e.first);
		std::sort((e.second)->begin(), (e.second)->end());
	} // for

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
	} // for

	return dist;
} // constructDistributions()

std::vector<rcptr<Graph>> GraphBuilder::constructClusters(
		const emdw::RVIds& vars, 
		std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses,
		std::map<emdw::RVIdType, rcptr<Factor>>& dist) const {
	
	std::map<emdw::RVIds, rcptr<Node>> nodes;
	std::vector<emdw::RVIds> subGraphScope(vars.size());

	// Step 1:
	for (unsigned i = 0; i < vars.size(); i++) {
		rcptr<DASS> aiDom = assocHypotheses[vars[i]];
		subGraphScope[i].push_back(vars[i]);

		for (unsigned j = i+1; j < vars.size(); j++) {
			rcptr<DASS> ajDom = assocHypotheses[vars[j]];

			DASS domIntersection;
			std::set_intersection( aiDom->begin(), aiDom->end(),
					ajDom->begin(), ajDom->end(),
					std::back_inserter(domIntersection));

			if (domIntersection.size() > 1) {
				emdw::RVIds scope = emdw::RVIds{vars[i], vars[j]};
				std::sort(scope.begin(), scope.end());
				
				if (!nodes[scope]) {
					rcptr<Factor> product = dist[vars[i]]->absorb(dist[vars[j]]);
					rcptr<DT> cancel = std::dynamic_pointer_cast<DT>(product);
					
					for (unsigned k = 1; k < domIntersection.size(); k++) {
						cancel->setEntry(scope, emdw::RVVals{ domIntersection[k], domIntersection[k]}, 0);
					} // for

					nodes[scope] = uniqptr<Node>(new Node(product));
				} // if
				subGraphScope[i].push_back(vars[j]);
			} // if
			domIntersection.clear();
		} // for
	} // for

	// Step 2:
	for (unsigned i = 0; i < vars.size(); i++) {
		for (unsigned j = 0; j < vars.size(); j++) {
			emdw::RVIds scopeIntersection;
			std::set_intersection( subGraphScope[i].begin(), subGraphScope[i].end(),
					subGraphScope[j].begin(), subGraphScope[j].end(),
					std::back_inserter(scopeIntersection));

			if ( scopeIntersection.size() > 0 ) {
				emdw::RVIds scopeUnion;
				std::set_union( subGraphScope[i].begin(), subGraphScope[i].end(),
						subGraphScope[j].begin(), subGraphScope[j].end(),
						std::back_inserter(scopeUnion));

				subGraphScope[i] = std::move(scopeUnion);
				scopeUnion.clear();
			} // if
			scopeIntersection.clear();
		} // for
	} // for

	// Step 3:
	std::sort(subGraphScope.begin(), subGraphScope.end());
	subGraphScope.erase( unique(subGraphScope.begin(), subGraphScope.end()), subGraphScope.end() );

	// Step 4:
	std::vector<std::vector<rcptr<Factor>>> subGraphs(subGraphScope.size());
	std::vector<rcptr<Graph>> graphs(subGraphScope.size());
	for (unsigned i = 0; i < subGraphScope.size(); i++) {
		for (auto const& n : nodes) {
			emdw::RVIds scopeIntersection;
			emdw::RVIds scope = f.first;

			std::set_intersection( subGraphScope[i].begin(), subGraphScope[i].end(),
					scope.begin(), scope.end(),
					std::back_inserter(scopeIntersection));

			if (scopeIntersection.size() > 1) subGraphs[i].push_back( factors[scope] );
		} // for

		if (subGraphs[i].size() > 1) {
			graphs[i] = uniqptr<Graph>(new Graph());

			for (auto& j : subGraphScope[i]) {
				for (unsigned k = 0; k < subGraphs[i].size(); k++) {
					
				} // for
			{ // for

		} // for
	} // for



	// Step 4:

	return graphs;
} // constructClusters()

