/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for a Simplified Conditonal Gaussian (Mixture) implementation.
 *************************************************************************/
#include <vector>
#include <iostream>
#include "sortindices.hpp"
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "matops.hpp"
#include "vecset.hpp"
#include "conditional_gaussian.hpp"

// Default operators
rcptr<FactorOperator> defaultInplaceNormalizerCG = uniqptr<FactorOperator>(new InplaceNormalizeCG());
rcptr<FactorOperator> defaultNormalizerCG = uniqptr<FactorOperator>(new NormalizeCG());
rcptr<FactorOperator> defaultInplaceAbsorberCG = uniqptr<FactorOperator>(new InplaceAbsorbCG());
rcptr<FactorOperator> defaultAbsorberCG = uniqptr<FactorOperator>(new AbsorbCG());
rcptr<FactorOperator> defaultInplaceCancellerCG = uniqptr<FactorOperator>(new InplaceCancelCG());
rcptr<FactorOperator> defaultCancellerCG = uniqptr<FactorOperator>(new CancelCG());
rcptr<FactorOperator> defaultMarginalizerCG = uniqptr<FactorOperator>(new MarginalizeCG());
rcptr<FactorOperator> defaultObserveReducerCG = uniqptr<FactorOperator>(new ObserveAndReduceCG());
rcptr<FactorOperator> defaultInplaceWeakDamperCG = uniqptr<FactorOperator>(new InplaceWeakDampingCG());

ConditionalGaussian::ConditionalGaussian(
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
	if (!inplaceNormalizer_) inplaceNormalizer_ = defaultInplaceNormalizerCG;
	if (!normalizer_) normalizer_ = defaultNormalizerCG;
	if (!inplaceAbsorber_) inplaceAbsorber_ = defaultInplaceAbsorberCG;
	if (!absorber_) absorber_ = defaultAbsorberCG;
	if (!inplaceCanceller_) inplaceCanceller_ = defaultInplaceCancellerCG;
	if (!canceller_) canceller_ = defaultCancellerCG;
	if (!marginalizer_) marginalizer_ = defaultMarginalizerCG;
	if (!observeAndReducer_) observeAndReducer_ = defaultObserveReducerCG;
	if (!inplaceDamper_) inplaceDamper_ = defaultInplaceWeakDamperCG;
} // Default Constructor

ConditionalGaussian::ConditionalGaussian(
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
	if (!inplaceNormalizer_) inplaceNormalizer_ = defaultInplaceNormalizerCG;
	if (!normalizer_) normalizer_ = defaultNormalizerCG;
	if (!inplaceAbsorber_) inplaceAbsorber_ = defaultInplaceAbsorberCG;
	if (!absorber_) absorber_ = defaultAbsorberCG;
	if (!inplaceCanceller_) inplaceCanceller_ = defaultInplaceCancellerCG;
	if (!canceller_) canceller_ = defaultCancellerCG;
	if (!marginalizer_) marginalizer_ = defaultMarginalizerCG;
	if (!observeAndReducer_) observeAndReducer_ = defaultObserveReducerCG;
	if (!inplaceDamper_) inplaceDamper_ = defaultInplaceWeakDamperCG;

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

ConditionalGaussian::~ConditionalGaussian() {} // Default Destructor

unsigned ConditionalGaussian::configure(unsigned) {
	std::cout << "NIY" << std::endl;
	return true;
} // configure()

unsigned ConditionalGaussian::classSpecificConfigure(
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
	this->~ConditionalGaussian();

	// .. and begin anew!
	new(this) ConditionalGaussian(
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
inline void ConditionalGaussian::inplaceNormalize(FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this);
	else dynamicInplaceApply(inplaceNormalizer_.get(), this);
} // inplaceNormalize()

inline uniqptr<Factor> ConditionalGaussian::normalize(FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor>(dynamicApply(procPtr, this));
	else return uniqptr<Factor>(dynamicApply(normalizer_.get(), this));
} // normalize()

//------------------Family 2: Absorbtion, Cancellation

inline void ConditionalGaussian::inplaceAbsorb(const Factor* rhsPtr, FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
	else dynamicInplaceApply(inplaceAbsorber_.get(), this, rhsPtr);
} // inplaceAbsorb()

inline uniqptr<Factor> ConditionalGaussian::absorb(const Factor* rhsPtr, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, rhsPtr));
	else return uniqptr<Factor> (dynamicApply(absorber_.get(), this, rhsPtr));
} // absorb()

inline void ConditionalGaussian::inplaceCancel(const Factor* rhsPtr, FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
	else dynamicInplaceApply(inplaceCanceller_.get(), this, rhsPtr);
} // inplaceCancel()

inline uniqptr<Factor> ConditionalGaussian::cancel(const Factor* rhsPtr, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, rhsPtr));
	else return uniqptr<Factor> (dynamicApply(canceller_.get(), this, rhsPtr));
} // cancel()

//------------------Family 4: Marginalization

inline uniqptr<Factor> ConditionalGaussian::marginalize(const emdw::RVIds& variablesToKeep, 
		bool presorted, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, variablesToKeep, presorted));
	else return uniqptr<Factor> (dynamicApply(marginalizer_.get(), this, variablesToKeep, presorted));
} // marginalize()

//------------------Family 4: ObserveAndReduce

inline uniqptr<Factor> ConditionalGaussian::observeAndReduce( const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, variables, assignedVals, presorted));
	else return uniqptr<Factor> (dynamicApply(observeAndReducer_.get(), this, variables, assignedVals, presorted));
} // observeAndReduce()


//------------------Family 4: Inplace Weak Damping

// TODO: Complete this!!!
double ConditionalGaussian::inplaceDampen(const Factor* oldMsg, double df, FactorOperator* procPtr) {
	if (procPtr) return dynamicInplaceApply(procPtr, this, oldMsg, df);
	else return dynamicInplaceApply(inplaceDamper_.get(), this, oldMsg, df); 
} // inplaceDampen()

//------------------Other required virtual methods

ConditionalGaussian* ConditionalGaussian::copy(const emdw::RVIds& newVars, bool presorted) const {

	if (newVars.size()) {
		rcptr<Factor> discrete = uniqptr<Factor>(discreteRV_->copy());

		std::map<unsigned, rcptr<Factor>> map; map.clear();

		for (auto& i : conditionalList_) {
			map[i.first] = uniqptr<Factor>((i.second)->copy());
		}

		return new ConditionalGaussian(discrete, 
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
	} // if
	
	return new ConditionalGaussian(*this);
} // copy()

ConditionalGaussian* ConditionalGaussian::vacuousCopy(const emdw::RVIds& selectedVars, bool presorted) const {
	return new ConditionalGaussian();
} // vacuousCopy()

bool ConditionalGaussian::isEqual(const Factor* rhsPtr) const { return true; } // isEqual()

unsigned ConditionalGaussian::noOfVars() const { return vars_.size(); } // noOfVars()

emdw::RVIds ConditionalGaussian::getVars() const { return vars_; } // getVars()

emdw::RVIdType ConditionalGaussian::getVar(unsigned varNo) const { return vars_[varNo]; } // getVar()

rcptr<Factor> ConditionalGaussian::getDiscretePrior() const { return discreteRV_; } // getDiscretePrior()

std::map<unsigned, rcptr<Factor>> ConditionalGaussian::getConditionalList() const { return conditionalList_; } // getConditionalList()

emdw::RVIds ConditionalGaussian::getContinuousVars() const {
	return (conditionalList_.begin()->second)->getVars();
} // getContinuousVars()

//TODO: Complete this!!!
std::istream& ConditionalGaussian::txtRead(std::istream& file) { return file; } // txtRead()

//TODO: Complete this!!
std::ostream& ConditionalGaussian::txtWrite(std::ostream& file) const { return file; } // txtWrite()

//==================================================FactorOperators======================================

//------------------Family 1: Normalization
//
const std::string& InplaceNormalizeCG::isA() const {
	static const std::string CLASSNAME("InplaceNormalizeCG");
	return CLASSNAME;
} // isA()

void InplaceNormalizeCG::inplaceProcess(ConditionalGaussian* lhsPtr) {
	ConditionalGaussian& lhs(*lhsPtr);

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

const std::string& NormalizeCG::isA() const {
	static const std::string CLASSNAME("NormalizeCG");
	return CLASSNAME;
} // isA()

Factor* NormalizeCG::process(const ConditionalGaussian* lhsPtr) {
	ConditionalGaussian* fPtr = new ConditionalGaussian(*lhsPtr);
	InplaceNormalizeCG ipNorm;
	
	try { 
		ipNorm.inplaceProcess(fPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()


//------------------Family 2: Absorption, Cancellation

const std::string& InplaceAbsorbCG::isA() const {
	static const std::string CLASSNAME("InplaceAbsorbCG");
	return CLASSNAME;
} // isA()

void InplaceAbsorbCG::inplaceProcess(ConditionalGaussian* lhsPtr, const Factor* rhsFPtr) {
	ConditionalGaussian& lhs(*lhsPtr);

	// Temporary variables
	rcptr<Factor> discretePrior;
	std::map<unsigned, rcptr<Factor>> map;
	rcptr<Factor> rhs = uniqptr<Factor>( rhsFPtr->copy() );
	const ConditionalGaussian* downCast;

	// An endless amount of options
	if (dynamic_cast<const ConditionalGaussian*>(rhsFPtr)) {
		downCast = dynamic_cast<const ConditionalGaussian*>(rhsFPtr);
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

const std::string& AbsorbCG::isA() const {
	static const std::string CLASSNAME("AbsorbCG");
	return CLASSNAME;
} // isA()

Factor* AbsorbCG::process(const ConditionalGaussian* lhsPtr, const Factor* rhsFPtr) {
	ConditionalGaussian* fPtr = new ConditionalGaussian(*lhsPtr);
	InplaceAbsorbCG ipAbsorb;
	
	try { 
		ipAbsorb.inplaceProcess(fPtr, rhsFPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()

const std::string& InplaceCancelCG::isA() const {
	static const std::string CLASSNAME("InplaceCancelCG");
	return CLASSNAME;
} // isA()

void InplaceCancelCG::inplaceProcess(ConditionalGaussian* lhsPtr, const Factor* rhsFPtr) {
	ConditionalGaussian& lhs(*lhsPtr);

	// Temporary variables
	rcptr<Factor> discretePrior;
	rcptr<Factor> mProj;
	std::map<unsigned, rcptr<Factor>> map;
	rcptr<Factor> rhs = uniqptr<Factor>( rhsFPtr->copy() );
	
	const ConditionalGaussian* downCast;
	const CanonicalGaussianMixture* gm;

	// An endless amount of options
	if (dynamic_cast<const ConditionalGaussian*>(rhsFPtr)) {
		downCast = dynamic_cast<const ConditionalGaussian*>(rhsFPtr);
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

const std::string& CancelCG::isA() const {
	static const std::string CLASSNAME("CancelCG");
	return CLASSNAME;
} // isA()


Factor* CancelCG::process(const ConditionalGaussian* lhsPtr, const Factor* rhsFPtr) {
	ConditionalGaussian* fPtr = new ConditionalGaussian(*lhsPtr);
	InplaceAbsorbCG ipCancel;
	
	try { 
		ipCancel.inplaceProcess(fPtr, rhsFPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()


//------------------Family 3: Marginalization

const std::string& MarginalizeCG::isA() const {
	static const std::string CLASSNAME("MarginalizeCG");
	return CLASSNAME;
} // isA()

Factor* MarginalizeCG::process(const ConditionalGaussian* lhsPtr, const emdw::RVIds& variablesToKeep, 
		bool presorted) {
	const ConditionalGaussian& lhs(*lhsPtr);

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
		map[i.first] = (i.second)->marginalize(continuousVars, presorted);
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

	return new ConditionalGaussian(discretePrior, 
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

const std::string& ObserveAndReduceCG::isA() const {
	static const std::string CLASSNAME("ObserveAndReduceCG");
	return CLASSNAME;
} // isA()

Factor* ObserveAndReduceCG::process(const ConditionalGaussian* lhsPtr, const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted) {
	const ConditionalGaussian& lhs(*lhsPtr);
	
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
		map[i.first]->inplaceNormalize();
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

	return new ConditionalGaussian(discretePrior, 
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

const std::string& InplaceWeakDampingCG::isA() const {
	static const std::string CLASSNAME("InplaceWeakDampingCG");
	return CLASSNAME;
} // isA()

// TODO: Complete this!!!
double InplaceWeakDampingCG::inplaceProcess(const ConditionalGaussian* lhsPtr, const Factor* rhsPtr, double df) {
	return 0.0;
} // inplaceProcess()
