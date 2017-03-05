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
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper) 
			: inplaceNormalizer_(inplaceNormalizer),
			normalizer_(normalizer),
			inplaceAbsorber_(inplaceAbsorber),
			absorber_(absorber),
			inplaceCanceller_(inplaceCanceller),
			canceller_(canceller),
			marginalizer_(marginalizer),
			observeAndReducer_(observerAndReducer),
			inplaceDamper_(inplaceDamper)
	{
	// Default operator intialisation
	if (!inplaceNormalizer_) inplaceNormalizer_ = defaultInplaceNormalizerLG;
	if (!normalizer_) normalizer_ = defaultNormalizerLG;
	if (!inplaceAbsorber_) inplaceAbsorber_ = defaultInplaceAbsorberLG;
	if (!absorber_) absorber_ = defaultAbsorberLG;
	if (!inplaceCanceller_) inplaceCanceller_ = defaultInplaceCancellerLG;
	if (!canceller_) canceller_ = defaultCancellerLG;
	if (!marginalizer_) marginalizer_ = defaultMarginalizerLG;
	if (!observeAndReducer_) observeAndReducer_ = defaultObserveReducerLG;
	if (!inplaceDamper_) inplaceDamper_ = defaultInplaceWeakDamperLG;
} // Default Constructor

LinearGaussian::LinearGaussian(
		const rcptr<Factor>& discreteRV,
		const std::map<unsigned, rcptr<Factor>>& conditionalList,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper) 
			: inplaceNormalizer_(inplaceNormalizer),
			normalizer_(normalizer),
			inplaceAbsorber_(inplaceAbsorber),
			absorber_(absorber),
			inplaceCanceller_(inplaceCanceller),
			canceller_(canceller),
			marginalizer_(marginalizer),
			observeAndReducer_(observerAndReducer),
			inplaceDamper_(inplaceDamper)
	{
	// Default operator intialisation
	if (!inplaceNormalizer_) inplaceNormalizer_ = defaultInplaceNormalizerLG;
	if (!normalizer_) normalizer_ = defaultNormalizerLG;
	if (!inplaceAbsorber_) inplaceAbsorber_ = defaultInplaceAbsorberLG;
	if (!absorber_) absorber_ = defaultAbsorberLG;
	if (!inplaceCanceller_) inplaceCanceller_ = defaultInplaceCancellerLG;
	if (!canceller_) canceller_ = defaultCancellerLG;
	if (!marginalizer_) marginalizer_ = defaultMarginalizerLG;
	if (!observeAndReducer_) observeAndReducer_ = defaultObserveReducerLG;
	if (!inplaceDamper_) inplaceDamper_ = defaultInplaceWeakDamperLG;

	// Get the variables
	emdw::RVIds vars;
	emdw::RVIds discreteVars = discreteRV->getVars(); 
	emdw::RVIds continuousVars = ( (conditionalList.begin())->second )->getVars();
	ASSERT( discreteVars.size() == 1, "discreteRV cannot be a distribution over " 
			<< discreteVars.size() << " variables, it must be non-vacuous over a single variable." );

	// Assign the discrete component
	vars.push_back(discreteVars[0]);
	isContinuous_[discreteVars[0]] = false;
	discreteRV_ = uniqptr<Factor>(discreteRV->copy());

	// Assign the continuous components
	for (auto& i : continuousVars) {
		vars.push_back(i);
		isContinuous_[i] = true;
	}

	for (auto& i : conditionalList) {
		ASSERT( continuousVars == (i.second)->getVars(), "All continuous distrubtions must be held the same variables " 
				<< continuousVars << " not" << (i.second)->getVars() );
		conditionalList_[i.first] = uniqptr<Factor>((i.second)->copy());
	}

	// Sort the variables
	std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
	vars_ = extract<unsigned>(vars, sorted);
} // Class Specific Constructor

LinearGaussian::~LinearGaussian() {} // Default Destructor

unsigned LinearGaussian::configure(unsigned) {
	std::cout << "NIY" << std::endl;
	return true;
} // configure()

unsigned LinearGaussian::classSpecificConfigure(
		const rcptr<Factor>& discreteRV,
		const std::map<unsigned, rcptr<Factor>>& conditionalList,
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
			discreteRV,
			conditionalList,
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
	else dynamicInplaceApply(inplaceCanceller_.get(), this, rhsPtr);
} // inplaceCancel()

inline uniqptr<Factor> LinearGaussian::cancel(const Factor* rhsPtr, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, rhsPtr));
	else return uniqptr<Factor> (dynamicApply(canceller_.get(), this, rhsPtr));
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
	else return uniqptr<Factor> (dynamicApply(observeAndReducer_.get(), this, variables, assignedVals, presorted));
} // observeAndReduce()


//------------------Family 4: Inplace Weak Damping

// TODO: Complete this!!!
double LinearGaussian::inplaceDampen(const Factor* oldMsg, double df, FactorOperator* procPtr) {
	if (procPtr) return dynamicInplaceApply(procPtr, this, oldMsg, df);
	else return dynamicInplaceApply(inplaceDamper_.get(), this, oldMsg, df); 
} // inplaceDampen()

//------------------Other required virtual methods

LinearGaussian* LinearGaussian::copy(const emdw::RVIds& newVars, bool presorted) const {
	rcptr<Factor> discrete = uniqptr<Factor>(discreteRV_->copy());

	std::map<unsigned, rcptr<Factor>> map; map.clear();

	for (auto& i : conditionalList_) {
		map[i.first] = uniqptr<Factor>((i.second)->copy());
	}

	return new LinearGaussian(discrete, 
				map,
				inplaceNormalizer_,
				normalizer_,
				inplaceAbsorber_,
				absorber_,
				inplaceCanceller_,
				canceller_,
				marginalizer_,
				observeAndReducer_,
				inplaceDamper_);
} // copy()

LinearGaussian* LinearGaussian::vacuousCopy(const emdw::RVIds& selectedVars, bool presorted) const {
	return new LinearGaussian();
} // vacuousCopy()

bool LinearGaussian::isEqual(const Factor* rhsPtr) const { return true; } // isEqual()

unsigned LinearGaussian::noOfVars() const { return vars_.size(); } // noOfVars()

emdw::RVIds LinearGaussian::getVars() const { return vars_; } // getVars()

emdw::RVIdType LinearGaussian::getVar(unsigned varNo) const { return vars_[varNo]; } // getVar()

rcptr<Factor> LinearGaussian::getDiscretePrior() const { return discreteRV_; } // getDiscretePrior()

std::map<unsigned, rcptr<Factor>> LinearGaussian::getConditionalList() const { return conditionalList_; } // getConditionalList()

emdw::RVIds LinearGaussian::getContinuousVars() const {
	return (conditionalList_.begin()->second)->getVars();
} // getContinuousVars()

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
	LinearGaussian& lhs(*lhsPtr);

	// Normalize the conditional Gaussians
	std::map<unsigned, rcptr<Factor>> map;
	for (auto& i : lhs.conditionalList_) map[i.first] = (i.second)->normalize();

	// Reconfigure the class
	lhs.classSpecificConfigure( (lhs.discreteRV_)->normalize(),
			        map,
				lhs.inplaceNormalizer_,
				lhs.normalizer_,
				lhs.inplaceAbsorber_,
				lhs.absorber_,
				lhs.inplaceCanceller_,
				lhs.canceller_,
				lhs.marginalizer_,
				lhs.observeAndReducer_,
				lhs.inplaceDamper_);
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
	LinearGaussian& lhs(*lhsPtr);

	// Temporary variables
	rcptr<Factor> discretePrior;
	std::map<unsigned, rcptr<Factor>> map;
	rcptr<Factor> rhs = uniqptr<Factor>( rhsFPtr->copy() );
	const LinearGaussian* downCast;

	// An endless amount of options
	if (dynamic_cast<const LinearGaussian*>(rhsFPtr)) {
		downCast = dynamic_cast<const LinearGaussian*>(rhsFPtr);
		ASSERT( (lhs.discreteRV_)->getVars() == (downCast->discreteRV_)->getVars(), 
				"The discrete components must have the same scope: "
				<< (lhs.discreteRV_)->getVars() << " != " << (downCast->discreteRV_)->getVars() );
		
		discretePrior = (lhs.discreteRV_)->absorb(rhs); // If the domains don't match everything should break here.
		for (auto& i : lhs.conditionalList_) map[i.first] = (i.second)->absorb( (downCast->conditionalList_)[i.first] );

	} else if (dynamic_cast<const GaussCanonical*>(rhsFPtr)) {
		discretePrior = uniqptr<Factor> ( (lhs.discreteRV_)->copy() );
		for (auto& i : lhs.conditionalList_) map[i.first] = uniqptr<Factor> ( (i.second)->absorb(rhs) );

	} else if (dynamic_cast<const CanonicalGaussianMixture*>(rhsFPtr)) {
		discretePrior = uniqptr<Factor> ( (lhs.discreteRV_)->copy() );
		for (auto& i : lhs.conditionalList_) map[i.first] = uniqptr<Factor> ( rhs->absorb(i.second) );

	} else if (dynamic_cast<const DiscreteTable<unsigned short>*>(rhsFPtr)) {
		ASSERT( (lhs.discreteRV_)->getVars() == rhs->getVars(), "The discrete distributions must have the same scope:" 
				<< (lhs.discreteRV_)->getVars() << " != " << rhs->getVars() );

		discretePrior = uniqptr<Factor> ( (lhs.discreteRV_)->absorb(rhs) );
		for (auto& i : lhs.conditionalList_) map[i.first] = uniqptr<Factor> ( (i.second)->copy() );
	}

	// Reconfigure the class
	lhs.classSpecificConfigure(discretePrior,
			        map,
				lhs.inplaceNormalizer_,
				lhs.normalizer_,
				lhs.inplaceAbsorber_,
				lhs.absorber_,
				lhs.inplaceCanceller_,
				lhs.canceller_,
				lhs.marginalizer_,
				lhs.observeAndReducer_,
				lhs.inplaceDamper_);
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
	LinearGaussian& lhs(*lhsPtr);

	// Temporary variables
	rcptr<Factor> discretePrior;
	rcptr<Factor> mProj;
	std::map<unsigned, rcptr<Factor>> map;
	rcptr<Factor> rhs = uniqptr<Factor>( rhsFPtr->copy() );
	
	const LinearGaussian* downCast;
	const CanonicalGaussianMixture* gm;

	// An endless amount of options
	if (dynamic_cast<const LinearGaussian*>(rhsFPtr)) {
		downCast = dynamic_cast<const LinearGaussian*>(rhsFPtr);
		ASSERT( (lhs.discreteRV_)->getVars() == (downCast->discreteRV_)->getVars(), 
				"The discrete components have the same single variable scope: " << (lhs.discreteRV_)->getVars() 
				<< " != " << (downCast->discreteRV_)->getVars() );
		
		discretePrior = (lhs.discreteRV_)->cancel(rhs); // If the domains don't match everything should break here.
		for (auto& i : lhs.conditionalList_) map[i.first] = (i.second)->cancel( (downCast->conditionalList_)[i.first] );

	} else if (dynamic_cast<const GaussCanonical*>(rhsFPtr)) {
		discretePrior = uniqptr<Factor> ( (lhs.discreteRV_)->copy() );
		for (auto& i : lhs.conditionalList_) map[i.first] = uniqptr<Factor> ( (i.second)->cancel(rhs) );

	} else if (dynamic_cast<const CanonicalGaussianMixture*>(rhsFPtr)) {
		gm = dynamic_cast<const CanonicalGaussianMixture*>(rhsFPtr);
		mProj = gm->momentMatch();

		discretePrior = uniqptr<Factor> ( (lhs.discreteRV_)->copy() );
		for (auto& i : lhs.conditionalList_) map[i.first] = uniqptr<Factor> ( (i.second)->cancel(mProj) );

	} else if (dynamic_cast<const DiscreteTable<unsigned short>*>(rhsFPtr)) {
		ASSERT( (lhs.discreteRV_)->getVars() == rhs->getVars(), 
			"The discrete components have the same single variable scope: "
			<< (lhs.discreteRV_)->getVars() << " != " << rhs->getVars() );

		discretePrior = uniqptr<Factor> ( (lhs.discreteRV_)->cancel(rhs) );
		for (auto& i : lhs.conditionalList_) map[i.first] = uniqptr<Factor> ( (i.second)->copy() );
	}

	// Reconfigure the class
	lhs.classSpecificConfigure(discretePrior,
			        map,
				lhs.inplaceNormalizer_,
				lhs.normalizer_,
				lhs.inplaceAbsorber_,
				lhs.absorber_,
				lhs.inplaceCanceller_,
				lhs.canceller_,
				lhs.marginalizer_,
				lhs.observeAndReducer_,
				lhs.inplaceDamper_);

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
	const LinearGaussian& lhs(*lhsPtr);

	//std::cout << "MarginalizeLG" << std::endl;

	// Temporary variables
	emdw::RVIds continuousVars; continuousVars.clear();
	emdw::RVIds discreteVar; discreteVar.clear();

	// Determine if the variables are discrete or not
	for (auto& i : variablesToKeep) {
		if ((lhs.isContinuous_)[i]) continuousVars.push_back(i);
		else discreteVar.push_back(i);
	} // for

	//std::cout << "continuousVars: " << continuousVars << std::endl;
	//std::cout << "discreteVar: " << discreteVar << std::endl;

	// Getting rid of continuous stuff usually happens
	rcptr<Factor> discretePrior = uniqptr<Factor>( (lhs.discreteRV_)->copy() );
	std::map<unsigned, rcptr<Factor>> map; map.clear();
	for (auto& i : lhs.conditionalList_) {
		//std::cout << "Before continuous marginalization: " << (i.second)->getVars() << std::endl;
		map[i.first] = (i.second)->marginalize(continuousVars, presorted)->normalize();
		//std::cout << "After continuous marginalization: " << map[i.first]->getVars() << std::endl;
	}

	// If you marginalize out all the continuous variables, you only have a discrete potential left.
	// TODO: Account for the mixtures' masses.
	if ( !(map.begin()->second)->noOfVars() && discreteVar.size() ) {
		return discretePrior->copy();
	} // if

	// If you're not keeping the discrete variable, you get a mixture.
	if (!discreteVar.size()) { 
		std::vector<rcptr<Factor>> mixtureComponents; mixtureComponents.clear();
		
		for(auto& i : map) {
			// Get the potential of the discrete variable
			rcptr<DiscreteTable<unsigned short>> dtConvert = 
				std::dynamic_pointer_cast<DiscreteTable<unsigned short>>(discretePrior);
			double potential = dtConvert->potentialAt(discretePrior->getVars(), 
					emdw::RVVals{ (unsigned short)(i.first) });

			// Get the factor
			rcptr<Factor> component = i.second;

			if (std::dynamic_pointer_cast<GaussCanonical>(component)) {
				// Adjust the mass
				std::dynamic_pointer_cast<GaussCanonical>(component)->adjustMass(potential);
				mixtureComponents.push_back(component);
			} else {
				rcptr<CanonicalGaussianMixture> cgmConvert = 
					std::dynamic_pointer_cast<CanonicalGaussianMixture>(component);
				std::vector<rcptr<Factor>> components = cgmConvert->getComponents();

				//std::cout << "CGM: " << components.size() << std::endl;

				// Adjust the mass
				for (rcptr<Factor> c : components) {
					std::dynamic_pointer_cast<GaussCanonical>(c)->adjustMass(potential);
					mixtureComponents.push_back(c);
				} // for
			}  // if
		} // for

		//std::cout << "mixturecomponents.size() = " << mixtureComponents.size() << std::endl;
		return new CanonicalGaussianMixture(variablesToKeep, mixtureComponents); // Default GM.
	} // if 

	return new LinearGaussian(discretePrior, 
			map,
			lhs.inplaceNormalizer_,
			lhs.normalizer_,
			lhs.inplaceAbsorber_,
			lhs.absorber_,
			lhs.inplaceCanceller_,
			lhs.canceller_,
			lhs.marginalizer_,
			lhs.observeAndReducer_,
			lhs.inplaceDamper_);
} // process()


//------------------Family 4: ObserveAndReduce

const std::string& ObserveAndReduceLG::isA() const {
	static const std::string CLASSNAME("ObserveAndReduceLG");
	return CLASSNAME;
} // isA()

Factor* ObserveAndReduceLG::process(const LinearGaussian* lhsPtr, const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted) {
	const LinearGaussian& lhs(*lhsPtr);
	
	// Temporary variables
	emdw::RVIds continuousVars, discreteVar;
	emdw::RVVals continuousVals, discreteVal;

	std::map<unsigned, rcptr<Factor>> map;
	rcptr<Factor> discretePrior;

	// Separate out continuous and discrete variables
	for (unsigned i = 0; i < variables.size(); i++) {
		if ((lhs.isContinuous_)[variables[i]]) {
			continuousVars.push_back(variables[i]);
			continuousVals.push_back(assignedVals[i]);
		} 
		else {
			discreteVar.push_back(variables[i]);
			discreteVal.push_back(assignedVals[i]);
		}
	}


	// Observing continuous things usually happens
	discretePrior = uniqptr<Factor> ( (lhs.discreteRV_)->copy() );
	for (auto& i : lhs.conditionalList_) {
		map[i.first] = (i.second)->observeAndReduce(continuousVars, continuousVals);
	}

	if ( discreteVar.size() ) {
		discretePrior = (lhs.discreteRV_)->observeAndReduce(discreteVar, discreteVal);
		rcptr<DiscreteTable<unsigned short>> dtConvert = 
			std::dynamic_pointer_cast<DiscreteTable<unsigned short>>(discretePrior);
		double potential = dtConvert->potentialAt(discreteVar, discreteVal);

		rcptr<GaussCanonical> gcConvert = std::dynamic_pointer_cast<GaussCanonical>( map[ (unsigned)(discreteVal[0]) ] );
		gcConvert->adjustMass(potential);

		return gcConvert->copy();
	}

	return new LinearGaussian(discretePrior, 
			map,
			lhs.inplaceNormalizer_,
			lhs.normalizer_,
			lhs.inplaceAbsorber_,
			lhs.absorber_,
			lhs.inplaceCanceller_,
			lhs.canceller_,
			lhs.marginalizer_,
			lhs.observeAndReducer_,
			lhs.inplaceDamper_);
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
