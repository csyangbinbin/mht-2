/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for the general cluster declared in Node.hpp
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
#include "node.hpp"

Node::Node(rcptr<Factor> cluster) : vars_(cluster->getVars()), prevFactor_(cluster->copy()) {
} // Constructor()

Node::~Node() {} // Default destructor()

void Node::addEdge(rcptr<Node> node, emdw::RVIds& sepset, rcptr<Factor> message) {
	sepsets_[node] = sepset;
	if (message) sentMsg_[node] = uniqptr<Factor> ( (factor_)->vacuousCopy(sepset, true) );
	else sentMsg_[node] = message;
} // addEdge()

emdw::RVIds Node::getVars() const {
	return vars_;
} // getVars()

rcptr<Factor> Node::getFactor() const {
	return uniqptr<Factor>( (factor_)->copy() );
} // getFactor()

std::vector<rcptr<Node>> Node::getAdjacentNodes() const {
	std::vector<rcptr<Node>> neighbours;
	for (auto it = sepsets_.begin(); it != sepsets_.end(); it++) neighbours.push_back(it->first);
	return neighbours;
} //getAdjacentNodes()

void Node::inplaceNormalize (FactorOperator* procPtr) {
	(factor_)->inplaceNormalize(procPtr);
} // inplaceNormalize()

uniqptr<Factor> Node::normalize (FactorOperator* procPtr) const {
	return (factor_)->normalize(procPtr);
} // normalize()

void Node::inplaceAbsorb (const Factor* rhsPtr, FactorOperator* procPtr) {
	(factor_)->inplaceAbsorb(rhsPtr, procPtr);
	vars_ = (factor_)->getVars();
} // inplaceAbsorb()

uniqptr<Factor> Node::absorb (const Factor* rhsPtr, FactorOperator* procPtr) const {
	return (factor_)->absorb(rhsPtr, procPtr);
} // absorb()

void Node::inplaceCancel (const Factor* rhsPtr, FactorOperator* procPtr) {
	(factor_)->inplaceCancel(rhsPtr, procPtr);
	vars_ = (factor_)->getVars();
} // inplaceCancel()

uniqptr<Factor> Node::cancel (const Factor* rhsPtr, FactorOperator* procPtr) const {
	return (factor_)->cancel(rhsPtr, procPtr);
} // cancel()

uniqptr<Factor> Node::marginalize(const emdw::RVIds& variablesToKeep, bool presorted, FactorOperator* procPtr) const {
	return (factor_)->marginalize(variablesToKeep, presorted, procPtr);
} // marginalize()

void Node::inplaceObserveAndReduce (const emdw::RVIds& variables, const emdw::RVVals& assignedVals, 
		bool presorted, FactorOperator* procPtr) {
	factor_ = (factor_)->observeAndReduce(variables, assignedVals, presorted, procPtr);
	vars_ = (factor_)->getVars();
} // inplaceObserveAndReduce()

uniqptr<Factor> Node::observeAndReduce (const emdw::RVIds& variables, const emdw::RVVals& assignedVals, 
		bool presorted, FactorOperator* procPtr) const {
	return (factor_)->observeAndReduce(variables, assignedVals, presorted, procPtr);
} // observeAndReduce()

template<typename T>
std::ostream& operator<<(std::ostream& file, const Node& n) {
	return (n.getFactor())->txtWrite(file);
} // operator<<
