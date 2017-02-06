/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for graph.hpp and node.hpp.
 *************************************************************************/
#include <iostream>
#include <map>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "anytype.hpp"
#include "emdw.hpp"
#include "oddsandsods.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "graph.hpp"
#include "node.hpp"
#include "transforms.hpp"
#include "utils.hpp"

using namespace emdw;
using namespace std;

class GraphTest : public testing::Test {
	protected:
		virtual void SetUp() {
		}

		virtual void TearDown() {
		}
};

TEST_F (GraphTest, GraphInitTest) {
	// Factors
	rcptr<Factor> gc_1 = uniqptr<Factor>(new GaussCanonical(emdw::RVIds{0, 1}));
	rcptr<Factor> gc_2 = uniqptr<Factor>(new GaussCanonical(emdw::RVIds{0, 2}));
	rcptr<Factor> gc_3 = uniqptr<Factor>(new GaussCanonical(emdw::RVIds{1, 2}));
	rcptr<Factor> gc_4 = uniqptr<Factor>(new GaussCanonical(emdw::RVIds{87, 54}));

	// Clusters
	rcptr<Node> psi_1 = uniqptr<Node>(new Node(gc_1));
	rcptr<Node> psi_2 = uniqptr<Node>(new Node(gc_2));
	rcptr<Node> psi_3 = uniqptr<Node>(new Node(gc_3));
	rcptr<Node> psi_4 = uniqptr<Node>(new Node(gc_4));

	// Graph
	rcptr<Graph> graph = uniqptr<Graph>(new Graph());
	graph->addEdge(psi_1, psi_2);
	graph->addEdge(psi_2, psi_3);
	graph->addEdge(psi_3, psi_1);

	// dfs
	std::map<rcptr<Node>, bool> marked = graph->depthFirstSearch();
}

TEST_F (GraphTest, NodeInitTest) {
	rcptr<Factor> gc = uniqptr<Factor>(new GaussCanonical());
	rcptr<Node> node = uniqptr<Node>(new Node(gc));
}
