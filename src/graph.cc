/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for the cluster graph declared in graph.hpp
 *************************************************************************/
#include <vector>
#include <map>
#include <iostream>
#include "sortindices.hpp"
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "matops.hpp"
#include "vecset.hpp"
#include "graph.hpp"
#include "node.hpp"

Graph::Graph () {} // Default constructor

Graph::Graph (const std::vector<rcptr<Node>>& nodes) {
	nodes_.clear(); 
	e_ = 0; n_ = 0;
	
	vars_.clear();
	contains_.clear();
	varMap_.clear();
	
	marked_.clear();
	converged_.clear();

	for (rcptr<Node> n : nodes) addNode(n);
} // Node vector constructor

Graph::~Graph () {} // Default destructor

void Graph::addNode (const rcptr<Node>& v) {
	emdw::RVIds tempUnion;
	emdw::RVIds vVars = v->getVars();

	// Append to the variable map.
	for(emdw::RVIdType i : vVars) {
		varMap_[i] = v;
		contains_[i] = true;
	} 
	
	// Add the node.
	nodes_.insert(v);
	n_ = nodes_.size();

	// Expand the scope.
	std::set_union(vVars.begin(), vVars.end(),
			vars_.begin(), vars_.end(),
			std::back_inserter(tempUnion));

	vars_ = std::move(tempUnion);
	tempUnion.clear(); vVars.clear();
} // addNode()

void Graph::addEdge (const rcptr<Node>& v, const rcptr<Node>& w) {
	emdw::RVIds l2i, r2i;
	emdw::RVIds sepset = sortedIntersection(v->getVars(), w->getVars(), l2i, r2i);

	ASSERT( sepset.size(), "v variables: " << v->getVars() << "w variables : " << w->getVars()
			<< " do not intersect, these clusters cannot share an edge." );
	
	addNode(v); addNode(w); // Hack, but easier than checking if nodes_ contains v, w

	v->addEdge(w, sepset);
	w->addEdge(v, sepset);
	e_++;
} // addEdge()

void Graph::depthFirstMessagePassing () {
	marked_.clear();

	for(rcptr<Node> v : nodes_) marked_[v] = false; 

	dfmp(*(nodes_.begin()));	

} // depthFirstMessagePassing()

void Graph::plotGraph () {
	marked_.clear();

	for(rcptr<Node> v : nodes_) marked_[v] = false; 

	dfp(*(nodes_.begin()));	
} // plotGraph()

void Graph::dfmp(const rcptr<Node>& v) {
	std::vector<std::weak_ptr<Node>> adjacent = v->getAdjacentNodes();

	marked_[v] = true;
	for (std::weak_ptr<Node> w : adjacent) bupReceiveMessage(v, w.lock());
	for (std::weak_ptr<Node> w : adjacent) if (!marked_[w]) dfmp(w.lock());

	v->cacheFactor(v->getFactor());
} // dfmp()

void Graph::bupReceiveMessage(const rcptr<Node>& v, const rcptr<Node>& w) {
	emdw::RVIds sepset = v->getSepset(w);

	// Divide the incoming message by the previous sent information
	rcptr<Factor> receivedMsg = w->getReceivedMessage(v);
	rcptr<Factor> incomingMsg = (w->marginalize(sepset))->cancel(receivedMsg);
	incomingMsg->inplaceNormalize();

	// Absorb and log the message
	v->inplaceAbsorb( incomingMsg.get() );
	v->inplaceNormalize();
	v->logMessage(w, incomingMsg);
} // bupReceiveMessage()

void Graph::dfp(const rcptr<Node>& v) {
	std::vector<std::weak_ptr<Node>> adjacent = v->getAdjacentNodes();
	marked_[v] = true;

	std::cout << "========================================" << std::endl;
	std::cout << *(v->getFactor()) << std::endl;
	std::cout << "----------------------------------------" << std::endl;
	std::cout << "Sepsets: " << std::endl;
	for (std::weak_ptr<Node> w: adjacent) {
		emdw::RVIds scope = (w.lock())->getVars();
		std::cout << v->getSepset(w.lock())[0] << " -> " << 
			"(" << scope[0] << "," << scope[1] << ")" << std::endl;
	}
	std::cout << "========================================" << std::endl;

	for (std::weak_ptr<Node> w : adjacent) if (!marked_[w]) dfp(w.lock());
} // dfp

std::vector<rcptr<Factor>> Graph::getFactors() const {
	std::vector<rcptr<Factor>> factors; factors.clear();

	for ( rcptr<Node> v : nodes_ ) factors.push_back( uniqptr<Factor>( (v->getFactor())->copy() ) );

	return factors;
} // getFactors()

bool Graph::containsVariable(const emdw::RVIdType i) const {
	return contains_[i];
} // containsVariable()

rcptr<Factor> Graph::getMarginalBelief(const emdw::RVIdType i) const {
	ASSERT( containsVariable(i), i << " is not in this graph!" );
	return ((varMap_[i]).lock())->marginalize(emdw::RVIds{i}, true);
} // getMarginalBelief()

unsigned Graph::getNoOfNodes() const {
	return n_;
} // getNoOfNodes()

unsigned Graph::getNoOfEdges() const {
	return e_;
} // getNoOfEdges()
