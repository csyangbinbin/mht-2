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

std::map<emdw::RVIdType, rcptr<Factor>> GraphBuilder::getMarginals(std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses) const {
	std::map<emdw::RVIdType, rcptr<Factor>> marginals;

	// Extract RVIds
	emdw::RVIds vars = extractRVIds(assocHypotheses);

	// Construct the distributions
	std::map<emdw::RVIdType, rcptr<Factor>> dist = constructDistributions(vars, assocHypotheses);

	// Construct the graphs
	marginals = constructClusters(vars, assocHypotheses, dist);

	return marginals;
} // getMarginals()

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
		
		sparseProbs[DASS{ (*aDom)[0] }] = 0.1;
		for (unsigned j = 1; j < aDom->size(); j++) sparseProbs[DASS{(*aDom)[j]}] = 1;
		
		dist[vars[i]] = uniqptr<Factor> (new DT(emdw::RVIds{vars[i]}, {aDom}, defProb_,
						sparseProbs, margin_, floor_, false,
						marginalizer_, inplaceNormalizer_, normalizer_ ) );

		sparseProbs.clear();
	} // for

	return dist;
} // constructDistributions()

std::map<emdw::RVIdType, rcptr<Factor>> GraphBuilder::constructClusters(
		const emdw::RVIds& vars, 
		std::map<emdw::RVIdType, rcptr<DASS>>& assocHypotheses,
		std::map<emdw::RVIdType, rcptr<Factor>>& dist) const {
	
	std::map<emdw::RVIdType, rcptr<Factor>> marginals; marginals.clear();
	std::vector<rcptr<Factor>> nodes; nodes.clear();
	std::vector<bool> connected(vars.size()); 
	
	// Step 1: Create the nodes and cancel conflicting hypotheses.
	for (unsigned i = 0; i < vars.size(); i++) {
		rcptr<DASS> aiDom = assocHypotheses[vars[i]];

		for (unsigned j = i+1; j < vars.size(); j++) {
			rcptr<DASS> ajDom = assocHypotheses[vars[j]];

			DASS domIntersection;
			std::set_intersection( aiDom->begin(), aiDom->end(),
					ajDom->begin(), ajDom->end(),
					std::back_inserter(domIntersection));

			if (domIntersection.size() > 1) {
				connected[i] = true; connected[j] = true;
				
				rcptr<Factor> product = dist[vars[i]]->absorb(dist[vars[j]]);
				rcptr<DT> cancel = std::dynamic_pointer_cast<DT>(product);
				emdw::RVIds scope = product->getVars();
				
				for (unsigned k = 1; k < domIntersection.size(); k++) {
					cancel->setEntry(scope, emdw::RVVals{ domIntersection[k], domIntersection[k]}, 0);
					product->inplaceNormalize();
				} // for
				nodes.push_back( product );
			} // if
			domIntersection.clear();
		} // for
	} // for

	// Step 2: Add in disjoint nodes
	for (unsigned i = 0; i < vars.size(); i++) {
		if (!connected[i]) {
			dist[vars[i]]->inplaceNormalize();
			nodes.push_back( uniqptr<Factor> (dist[vars[i]]->copy()) );
		} // if
	} // for

	// Step 3: Create the cluster graph
	rcptr<ClusterGraph> clusterGraph = uniqptr<ClusterGraph>(new ClusterGraph(nodes));
	std::map<Idx2, rcptr<Factor>> msgs; msgs.clear();
	MessageQueue msgQ; msgQ.clear();

	// Step 4: Pass messages until convergence
	unsigned nMsg = loopyBU_CG(*clusterGraph, msgs, msgQ, 0.0);

	// Step 5: Extract the marginals
	for (emdw::RVIdType i : vars)  marginals[i] = queryLBU_CG(*clusterGraph, msgs, emdw::RVIds{i} )->normalize();

	return marginals;
} // constructClusters()

