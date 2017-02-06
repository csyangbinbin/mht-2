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

Graph::Graph (std::vector<rcptr<Node>> nodes) {
	for (rcptr<Node> n : nodes) addNode(n);
} // Node vector constructor

Graph::~Graph () {} // Default destructor

void Graph::addNode (rcptr<Node> v) {
	nodes_.insert(v);
	n_ = nodes_.size();
} // addNode()

void Graph::addEdge (rcptr<Node> v, rcptr<Node> w) {
	emdw::RVIds l2i, r2i;
	emdw::RVIds sepset = sortedIntersection(v->getVars(), w->getVars(), l2i, r2i);

	ASSERT( sepset.size(), "v variables: " << v->getVars() << "w variables : " << w->getVars()
			<< " do not intersect, these clusters cannot share an edge." );
	
	addNode(v); addNode(w); // Hack, but easier than checking if nodes_ contains v, w

	v->addEdge(w, sepset);
	w->addEdge(v, sepset);
	e_++;
} // addEdge()

std::map<rcptr<Node>, bool> Graph::depthFirstSearch () {
	marked_.clear();

	for(rcptr<Node> v : nodes_) marked_[v] = false; 

	dfs(*(nodes_.begin()));
	
	return marked_;
} // depthFirstSearch()

void Graph::dfs(rcptr<Node> v) {
	std::cout << v->getVars() << std::endl;

	marked_[v] = true;

	for (rcptr<Node> w : v->getAdjacentNodes()) {
		if (!marked_[w]) dfs(w);
	}

} // dfs()

unsigned Graph::getNoOfNodes() const {
	return n_;
} // getNoOfNodes()

unsigned Graph::getNoOfEdges() const {
	return e_;
} // getNoOfEdges()
