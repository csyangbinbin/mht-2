/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for a Simplified Linear Gaussian implementation.
 *************************************************************************/
#include <vector>
#include <iostream>
#include "sortindices.hpp"
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "matops.hpp"
#include "vecset.hpp"
#include "linear_gaussian.hpp"

// Default operators
rcptr<FactorOperator> defaultInplaceNormalizerLG = uniqptr<FactorOperator>(new InplaceNormalizeLG());
rcptr<FactorOperator> defaultNormalizerLG = uniqptr<FactorOperator>(new NormalizeLG());
rcptr<FactorOperator> defaultInplaceAbsorberLG = uniqptr<FactorOperator>(new InplaceAbsorbLG());
rcptr<FactorOperator> defaultAbsorberLG = uniqptr<FactorOperator>(new AbsorbLG());
rcptr<FactorOperator> defaultInplaceCancellerLG = uniqptr<FactorOperator>(new InplaceCancelLG());
rcptr<FactorOperator> defaultCancellerLG = uniqptr<FactorOperator>(new CancelLG());
rcptr<FactorOperator> defaultMarginalizerLG = uniqptr<FactorOperator>(new MarginalizeLG());
rcptr<FactorOperator> defaultObserveReducerLG = uniqptr<FactorOperator>(new ObserveAndReduceLG());
rcptr<FactorOperator> defaultInplaceWeakDamperLG = uniqptr<FactorOperator>(new InplaceWeakDampingLG());

LinearGaussian::LinearGaussian(
		const emdw::RVIds& vars,
		bool presorted,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper)
		{
	
	// Default operator intialisation
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerLG;
	if (!normalizer) normalizer_ = defaultNormalizerLG;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorberLG;
	if (!absorber) absorber_ = defaultAbsorberLG;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCancellerLG;
	if (!canceller) canceller_ = defaultCancellerLG;
	if (!marginalizer) marginalizer_ = defaultMarginalizerLG;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducerLG;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamperLG;

	// Ensure the higher level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;
	} else {
		std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
		vars_ = extract<unsigned>(vars, sorted);
	}
} // Default Constructor

LinearGaussian::~LinearGaussian() {} // Default Destructor

unsigned LinearGaussian::classSpecificConfigure(
		const emdw::RVIds& vars,
		bool presorted,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {
	// Destroy existing ...
	this->~LinearGaussian();

	// .. and begin anew!
	new(this) LinearGaussian(
			vars,
			presorted,
			inplaceNormalizer,
			normalizer,
			inplaceAbsorber,
			absorber,
			inplaceCanceller,
			canceller,
			marginalizer,
			observerAndReducer,
			inplaceDamper);
	
	return 1;
} // classSpecificConfigure()


//------------------Family 1: Normalization

inline void LinearGaussian::inplaceNormalize(FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this);
	else dynamicInplaceApply(inplaceNormalizer_.get(), this);
} // inplaceNormalize()

inline uniqptr<Factor> LinearGaussian::normalize(FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor>(dynamicApply(procPtr, this));
	else return uniqptr<Factor>(dynamicApply(normalizer_.get(), this));
} // normalize()

//------------------Family 2: Absorbtion, Cancellation

inline void LinearGaussian::inplaceAbsorb(const Factor* rhsPtr, FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
	else dynamicInplaceApply(inplaceAbsorber_.get(), this, rhsPtr);
} // inplaceAbsorb()

inline uniqptr<Factor> LinearGaussian::absorb(const Factor* rhsPtr, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, rhsPtr));
	else return uniqptr<Factor> (dynamicApply(absorber_.get(), this, rhsPtr));
} // absorb()

inline void LinearGaussian::inplaceCancel(const Factor* rhsPtr, FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
	else dynamicInplaceApply(inplaceAbsorber_.get(), this, rhsPtr);
} // inplaceCancel()

inline uniqptr<Factor> LinearGaussian::cancel(const Factor* rhsPtr, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, rhsPtr));
	else return uniqptr<Factor> (dynamicApply(absorber_.get(), this, rhsPtr));
} // cancel()

//------------------Family 4: Marginalization

inline uniqptr<Factor> LinearGaussian::marginalize(const emdw::RVIds& variablesToKeep, 
		bool presorted, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, variablesToKeep, presorted));
	else return uniqptr<Factor> (dynamicApply(marginalizer_.get(), this, variablesToKeep, presorted));
} // marginalize()

//------------------Family 4: ObserveAndReduce

inline uniqptr<Factor> LinearGaussian::observeAndReduce( const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, variables, assignedVals, presorted));
	else return uniqptr<Factor> (dynamicApply(marginalizer_.get(), this, variables, assignedVals, presorted));
} // observeAndReduce()


//------------------Family 4: Inplace Weak Damping

// TODO: Complete this!!!
double LinearGaussian::inplaceDampen(const Factor* oldMsg, double df, FactorOperator* procPtr) {
	if (procPtr) return dynamicInplaceApply(procPtr, this, oldMsg, df);
	else return dynamicInplaceApply(inplaceDamper_.get(), this, oldMsg, df); 
} // inplaceDampen()

//------------------Other required virtual methods

LinearGaussian* LinearGaussian::copy(const emdw::RVIds& newVars, bool presorted) const {
	return new LinearGaussian(*this);
} // copy()

LinearGaussian* LinearGaussian::vacuousCopy(const emdw::RVIds& selectedVars, bool presorted) const {
	return new LinearGaussian();
} // vacuousCopy()

bool LinearGaussian::isEqual(const Factor* rhsPtr) const { return true; } // isEqual()

unsigned LinearGaussian::noOfVars() const { return vars_.size(); } // noOfVars()

emdw::RVIds LinearGaussian::getVars() const { return vars_; } // getVars()

emdw::RVIdType LinearGaussian::getVar(unsigned varNo) const { return vars_[varNo]; } // getVar()

//TODO: Complete this!!!
std::istream& LinearGaussian::txtRead(std::istream& file) { return file; } // txtRead()

//TODO: Complete this!!
std::ostream& LinearGaussian::txtWrite(std::ostream& file) const { return file; } // txtWrite()

//==================================================FactorOperators======================================

//------------------Family 1: Normalization
//
const std::string& InplaceNormalizeLG::isA() const {
	static const std::string CLASSNAME("InplaceNormalizeLG");
	return CLASSNAME;
} // isA()

void InplaceNormalizeLG::inplaceProcess(LinearGaussian* lhsPtr) {
} // inplaceProcess()

const std::string& NormalizeLG::isA() const {
	static const std::string CLASSNAME("NormalizeLG");
	return CLASSNAME;
} // isA()

Factor* NormalizeLG::process(const LinearGaussian* lhsPtr) {
	LinearGaussian* fPtr = new LinearGaussian(*lhsPtr);
	InplaceNormalizeLG ipNorm;
	
	try { 
		ipNorm.inplaceProcess(fPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()


//------------------Family 2: Absorption, Cancellation

const std::string& InplaceAbsorbLG::isA() const {
	static const std::string CLASSNAME("InplaceAbsorbLG");
	return CLASSNAME;
} // isA()

void InplaceAbsorbLG::inplaceProcess(LinearGaussian* lhsPtr, const Factor* rhsFPtr) {
} // inplaceProcess()

const std::string& AbsorbLG::isA() const {
	static const std::string CLASSNAME("AbsorbLG");
	return CLASSNAME;
} // isA()

Factor* AbsorbLG::process(const LinearGaussian* lhsPtr, const Factor* rhsFPtr) {
	LinearGaussian* fPtr = new LinearGaussian(*lhsPtr);
	InplaceAbsorbLG ipAbsorb;
	
	try { 
		ipAbsorb.inplaceProcess(fPtr, rhsFPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()

const std::string& InplaceCancelLG::isA() const {
	static const std::string CLASSNAME("InplaceCancelLG");
	return CLASSNAME;
} // isA()

void InplaceCancelLG::inplaceProcess(LinearGaussian* lhsPtr, const Factor* rhsFPtr) {
} // inplaceCancel()

const std::string& CancelLG::isA() const {
	static const std::string CLASSNAME("CancelLG");
	return CLASSNAME;
} // isA()


Factor* CancelLG::process(const LinearGaussian* lhsPtr, const Factor* rhsFPtr) {
	LinearGaussian* fPtr = new LinearGaussian(*lhsPtr);
	InplaceAbsorbLG ipCancel;
	
	try { 
		ipCancel.inplaceProcess(fPtr, rhsFPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()


//------------------Family 3: Marginalization

const std::string& MarginalizeLG::isA() const {
	static const std::string CLASSNAME("MarginalizeLG");
	return CLASSNAME;
} // isA()

Factor* MarginalizeLG::process(const LinearGaussian* lhsPtr, const emdw::RVIds& variablesToKeep, 
		bool presorted) {
	return new GaussCanonical();
} // process()


//------------------Family 4: ObserveAndReduce

const std::string& ObserveAndReduceCGM::isA() const {
	static const std::string CLASSNAME("ObserveAndReduceLG");
	return CLASSNAME;
} // isA()

Factor* ObserveAndReduceLG::process(const LinearGaussian* lhsPtr, const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted) {
	return new GaussCanonical();
} // process()


//------------------Family 5: Damping

const std::string& InplaceWeakDampingLG::isA() const {
	static const std::string CLASSNAME("InplaceWeakDampingLG");
	return CLASSNAME;
} // isA()

// TODO: Complete this!!!
double InplaceWeakDampingLG::inplaceProcess(const LinearGaussian* lhsPtr, const Factor* rhsPtr, double df) {
	return 0.0;
} // inplaceProcess()
