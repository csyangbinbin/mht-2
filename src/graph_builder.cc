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
rcptr<FactorOperator> defaultNormalizerGB = uniqptr<FactorOperator>(new DiscreteTable_Normalize<unsigned short>);
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
		
		sparseProbs[DASS{ (*aDom)[0] }] = 0.5;
		for (unsigned j = 0; j < aDom->size(); j++) sparseProbs[DASS{(*aDom)[j]}] = 1;
		
		dist[vars[i]] = uniqptr<Factor> (new DT(emdw::RVIds{vars[i]}, {aDom}, defProb_,
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
	std::vector<bool> connected(vars.size()); 
	
	// Step 1: Create the nodes and cancel conflicting hypotheses.
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
				connected[i] = true; connected[j] = true;
				emdw::RVIds scope = emdw::RVIds{vars[i], vars[j]};
				std::sort(scope.begin(), scope.end());
				
				if (!nodes[scope]) {
					rcptr<Factor> product = dist[vars[i]]->absorb(dist[vars[j]]);
					rcptr<DT> cancel = std::dynamic_pointer_cast<DT>(product);
					
					for (unsigned k = 1; k < domIntersection.size(); k++) {
						cancel->setEntry(scope, emdw::RVVals{ domIntersection[k], domIntersection[k]}, 0);
						product->inplaceNormalize();
					} // for

					nodes[scope] = uniqptr<Node>(new Node(product));
				} // if
				subGraphScope[i].push_back(vars[j]);
			} // if
			domIntersection.clear();
		} // for
	} // for

	// Step 2: Add in disjoint nodes
	for (unsigned i = 0; i < vars.size(); i++) {
		if (!connected[i]) {
			dist[vars[i]]->inplaceNormalize();
			nodes[emdw::RVIds{vars[i]}] = uniqptr<Node>(new Node(dist[vars[i]]));
		} // if
	} // for

	// Step 3: Determine all disjoint association networks
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

	// Step 4: Determine disjoint subgraphs
	std::sort(subGraphScope.begin(), subGraphScope.end());
	subGraphScope.erase( unique(subGraphScope.begin(), subGraphScope.end()), subGraphScope.end() );

	// Step 5: Link up subgraph nodes
	std::vector<std::vector<rcptr<Node>>> subGraph(subGraphScope.size());
	std::vector<rcptr<Graph>> graph(subGraphScope.size());
	for (unsigned i = 0; i < subGraphScope.size(); i++) {
		// Step 5.1: Determine whether a node is in a subgraph.
		for (auto const& n : nodes) {
			emdw::RVIds scopeIntersection;
			emdw::RVIds scope = n.first;

			std::set_intersection( subGraphScope[i].begin(), subGraphScope[i].end(),
					scope.begin(), scope.end(),
					std::back_inserter(scopeIntersection));

			if (scopeIntersection.size() > 0) subGraph[i].push_back( nodes[scope] );
			scopeIntersection.clear();
		} // for

		// Step 5.2 : Link up the nodes in a chain.
		if (subGraph[i].size() > 1) {
			graph[i] = uniqptr<Graph>(new Graph());
			for (emdw::RVIdType j : subGraphScope[i]) {
				int prev = -1;
				for (unsigned k = 0; k < subGraph[i].size(); k++) {
					bool intersects = false;
					emdw::RVIds nodeVars = subGraph[i][k]->getVars();

					if (j == nodeVars[0] || j == nodeVars[1]) intersects = true;
					if (intersects && prev != -1) graph[i]->addEdge(subGraph[i][prev], subGraph[i][k]);
					if (intersects) prev = (int) k;

					nodeVars.clear();
				} // for
			} // for
		} else graph[i] = uniqptr<Graph>(new Graph( { subGraph[i][0] } )); // If there is only a single node.
	} // for

	return graph;
} // constructClusters()

