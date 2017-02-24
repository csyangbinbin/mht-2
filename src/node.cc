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

Node::Node(const rcptr<Factor>& factor) {
	factor_ = uniqptr<Factor>( factor->copy() );
	vars_ = factor_->getVars();
	
	adjacent_.clear();
	sepsets_.clear();

	prevFactor_ = uniqptr<Factor>( factor->copy() );
	recMsg_.clear();
} // Constructor()

Node::~Node() {
} // Default destructor()

void Node::addEdge(const rcptr<Node>& w, const emdw::RVIds& sepset, const rcptr<Factor>& message) {
	adjacent_.push_back(w);
	sepsets_.push_back(sepset);
} // addEdge()

void Node::logMessage(const rcptr<Node>& w, const rcptr<Factor>& message) {
	//recMsg_[w] = message;
} // logMessage()

void Node::cacheFactor(const rcptr<Factor>& factor) {
	prevFactor_ = uniqptr<Factor>(factor->copy());
} // cacheFactor()

emdw::RVIds Node::getVars() const {
	return vars_;
} // getVars()

rcptr<Factor> Node::getFactor() const {
	return uniqptr<Factor>( (factor_)->copy() );
} // getFactor()

rcptr<Factor> Node::getCachedFactor() const {
	return uniqptr<Factor> ( (prevFactor_)->copy() );
} // getCachedFactor()

emdw::RVIds Node::getSepset(const rcptr<Node>& w) {
	//return sepsets_[w];
	return emdw::RVIds{0};
} // getSepset()

rcptr<Factor> Node::getReceivedMessage(const rcptr<Node>& w) {
	//return recMsg_[w];
	return uniqptr<Factor>(factor_->copy());
} // getSentMessage()

std::vector<rcptr<Node>> Node::getAdjacentNodes() const {
	return std::vector<rcptr<Node>>(2);
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
