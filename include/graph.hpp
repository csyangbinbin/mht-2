/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for a graph.
 *************************************************************************/
#ifndef GRAPH_HPP
#define GRAPH_HPP

#include "factor.hpp"
#include "factoroperator.hpp"
#include "emdw.hpp"
#include "anytype.hpp"
#include "node.hpp"

class Graph;

/**
 * @brief A handcrafted cluster graph.
 *
 * @author SCJ Robertson
 * @since 06/02/17
 */
class Graph {

	public:
		/**
		 * @brief Default constructor.
		 */
		Graph();

		/**
		 * @brief Node vector constructor.
		 *
		 * @param nodes A collection of nodes.
		 */
		Graph(std::vector<rcptr<Node>> nodes);

		/**
		 * @brief Default destructor.
		 */
		~Graph();

	public:
		/**
		 * @brief Add a new cluster to the graph
		 *
		 * @param v A node within the graph.
		 */
		void addNode (rcptr<Node> v);

		/**
		 * @brief Add an edge between two clusters
		 *
		 * @param v A node within the graph, to be
		 * connected with w.
		 *
		 * @param w A node within the graph, to 
		 * be connected to v.
		 */
		void addEdge (rcptr<Node> v, rcptr<Node> w);

	public:
		/**
		 * @brief Depth First Search
		 */
		std::map<rcptr<Node>, bool> depthFirstSearch();

	private:
		/**
		 * @brief Recursive Depth Frist search
		 */
		void dfs(rcptr<Node> v);
	
	public:
		/**
		 * @brief Get number of nodes in the graph.
		 */
		unsigned getNoOfNodes() const;

		/**
		 * @brief Get the number of edges in the graph.
		 */
		unsigned getNoOfEdges() const;

	private:
		std::set<rcptr<Node>> nodes_;
		mutable unsigned n_, e_;

		std::map<rcptr<Node>, bool> marked_;

}; // Graph

#endif // GRAPH_HPP
