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
		Graph(const std::vector<rcptr<Node>>& nodes);

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
		void addNode (const rcptr<Node>& v);

		/**
		 * @brief Add an edge between two clusters
		 *
		 * The sepsets on this edge will be an exact intersection
		 * of the adjacent clusters' scopes, because I'm so
		 * lazy.
		 *
		 * @param v A node within the graph, to be
		 * connected with w.
		 *
		 * @param w A node within the graph, to 
		 * be connected to v.
		 */
		void addEdge (const rcptr<Node>& v, const rcptr<Node>& w);

	public:
		/**
		 * @brief Depth First Message Passing
		 *
		 * Depth First Search, a naive attempt at a 
		 * message passing schedule which at least ensure every node
		 * in the graph is visited.
		 */
		void depthFirstMessagePassing();

	private:
		/**
		 * @brief Recursive Depth First Message Passing
		 *
		 * @param v The current node in the graph.
		 */
		void dfmp(const rcptr<Node>& v);

		/**
		 * @brief Incoporate a message into a given cluster.
		 *
		 * Incoporate a message into a given cluster under the BUP 
		 * framework.
		 *
		 * @param v The receiving node.
		 *
		 * @param w The transimitting node.
		 */
		void bupReceiveMessage(const rcptr<Node>& v, const rcptr<Node>& w);

		/**
		 * @brief Recursive Depth First Search Plot, as called by plotGraph.
		 *
		 * @param v The current node in the graph.
		 */
		void dfp(const rcptr<Node>& v);

	public:
		/**
		 * @brief Returns the factors in the graph.
		 */
		std::vector<rcptr<Factor>> getFactors() const;

		/**
		 * @brief Return number of nodes in the graph.
		 */
		unsigned getNoOfNodes() const;

		/**
		 * @brief Return the number of edges in the graph.
		 */
		unsigned getNoOfEdges() const;

	public:
		/**
		 * @brief Print the Graph in DFS order.
		 */
		void plotGraph();

	private:
		std::set<rcptr<Node>> nodes_;
		mutable unsigned n_, e_;
		emdw::RVIds vars_;

		std::map<std::weak_ptr<Node>, bool, std::owner_less<std::weak_ptr<Node>>> marked_;
		std::map<std::weak_ptr<Node>, bool, std::owner_less<std::weak_ptr<Node>>> converged_;

}; // Graph

#endif // GRAPH_HPP
