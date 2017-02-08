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

Graph::Graph (const std::vector<rcptr<Node>> nodes) {
	for (rcptr<Node> n : nodes) addNode(n);
} // Node vector constructor

Graph::~Graph () {} // Default destructor

void Graph::addNode (const rcptr<Node> v) {
	nodes_.insert(v);
	n_ = nodes_.size();
} // addNode()

void Graph::addEdge (const rcptr<Node> v, const rcptr<Node> w) {
	emdw::RVIds l2i, r2i;
	emdw::RVIds sepset = sortedIntersection(v->getVars(), w->getVars(), l2i, r2i);

	ASSERT( sepset.size(), "v variables: " << v->getVars() << "w variables : " << w->getVars()
			<< " do not intersect, these clusters cannot share an edge." );
	
	addNode(v); addNode(w); // Hack, but easier than checking if nodes_ contains v, w

	v->addEdge(w, sepset);
	w->addEdge(v, sepset);
	e_++;
} // addEdge()

void Graph::depthFirstSearch () {
	marked_.clear();

	for(rcptr<Node> v : nodes_) marked_[v] = false; 

	dfs(*(nodes_.begin()));	

} // depthFirstSearch()

void Graph::dfs(const rcptr<Node> v) {

	std::cout << v->getVars() << std::endl;

	std::vector<rcptr<Node>> adjacent = v->getAdjacentNodes();

	marked_[v] = true;
	//for (rcptr<Node> w : adjacent) bupReceiveMessage(v, w);
	for (rcptr<Node> w : adjacent) if (!marked_[w]) dfs(w);

	v->cacheFactor(v->getFactor());
} // dfs()

void Graph::bupReceiveMessage(const rcptr<Node> v, const rcptr<Node> w) {
	emdw::RVIds sepset = v->getSepset(w);

	// Divide the incoming message by the previous sent information
	rcptr<Factor> receivedMsg = w->getReceivedMessage(v);
	rcptr<Factor> incomingMsg = (w->marginalize(sepset))->cancel(receivedMsg);

	// Absorb and log the message
	v->inplaceAbsorb(incomingMsg->copy());
	v->inplaceNormalize();
	v->logMessage(w, incomingMsg);
} // bupReceiveMessage()

unsigned Graph::getNoOfNodes() const {
	return n_;
} // getNoOfNodes()

unsigned Graph::getNoOfEdges() const {
	return e_;
} // getNoOfEdges()
