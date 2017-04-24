/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for the Canonical Gaussian Mixture implementation.
 *************************************************************************/
#include <vector>
#include <map>
#include <iostream>
#include <math.h>
#include "sortindices.hpp"
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "matops.hpp"
#include "vecset.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"

// Default operators
rcptr<FactorOperator> defaultInplaceNormalizerCGM = uniqptr<FactorOperator>(new InplaceNormalizeCGM());
rcptr<FactorOperator> defaultNormalizerCGM = uniqptr<FactorOperator>(new NormalizeCGM());
rcptr<FactorOperator> defaultInplaceAbsorberCGM = uniqptr<FactorOperator>(new InplaceAbsorbCGM());
rcptr<FactorOperator> defaultAbsorberCGM = uniqptr<FactorOperator>(new AbsorbCGM());
rcptr<FactorOperator> defaultInplaceCancellerCGM = uniqptr<FactorOperator>(new InplaceCancelCGM());
rcptr<FactorOperator> defaultCancellerCGM = uniqptr<FactorOperator>(new CancelCGM());
rcptr<FactorOperator> defaultMarginalizerCGM = uniqptr<FactorOperator>(new MarginalizeCGM());
rcptr<FactorOperator> defaultObserveReducerCGM = uniqptr<FactorOperator>(new ObserveAndReduceCGM());
rcptr<FactorOperator> defaultInplaceWeakDamperCGM = uniqptr<FactorOperator>(new InplaceWeakDampingCGM());

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const emdw::RVIds& vars,
		bool presorted,
		const unsigned maxComponents,
		const double threshold,
		const double unionDistance,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper)
			: vars_(vars.size()),
			maxComp_(maxComponents),
			threshold_(threshold),
			unionDistance_(unionDistance),
			inplaceNormalizer_(inplaceNormalizer),
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
	if (!inplaceNormalizer_) { inplaceNormalizer_ = defaultInplaceNormalizerCGM; }
	if (!normalizer_) { normalizer_ = defaultNormalizerCGM; }
	if (!inplaceAbsorber_) { inplaceAbsorber_ = defaultInplaceAbsorberCGM; }
	if (!absorber_) { absorber_ = defaultAbsorberCGM; }
	if (!inplaceCanceller_) { inplaceCanceller_ = defaultInplaceCancellerCGM; }
	if (!canceller_) { canceller_ = defaultCancellerCGM; }
	if (!marginalizer_) { marginalizer_ = defaultMarginalizerCGM; }
	if (!observeAndReducer_) { observeAndReducer_ = defaultObserveReducerCGM; }
	if (!inplaceDamper_) { inplaceDamper_ = defaultInplaceWeakDamperCGM; }

	// Ensure the higher level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;
	} else {
		std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
		vars_ = extract<unsigned>(vars, sorted);
	}

	// Create a mixture with a single vacuous component.
	comps_.clear();
	comps_.push_back( uniqptr<Factor> ( new GaussCanonical(vars_, true ) ) );
	N_ = 1;
} // Default Constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const emdw::RVIds& vars,
		const std::vector<double>& weights,
		const std::vector<ColVector<double>>& means,
		const std::vector<Matrix<double>>& covs,
		bool presorted,
		const unsigned maxComponents,
		const double threshold,
		const double unionDistance,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper)
			: vars_(vars.size()),
			comps_(weights.size()),
			N_(weights.size()),
			maxComp_(maxComponents),
			threshold_(threshold),
			unionDistance_(unionDistance),
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
	if (!inplaceNormalizer_) { inplaceNormalizer_ = defaultInplaceNormalizerCGM; }
	if (!normalizer_) { normalizer_ = defaultNormalizerCGM; }
	if (!inplaceAbsorber_) { inplaceAbsorber_ = defaultInplaceAbsorberCGM; }
	if (!absorber_) { absorber_ = defaultAbsorberCGM; }
	if (!inplaceCanceller_) { inplaceCanceller_ = defaultInplaceCancellerCGM; }
	if (!canceller_) { canceller_ = defaultCancellerCGM; }
	if (!marginalizer_) { marginalizer_ = defaultMarginalizerCGM; }
	if (!observeAndReducer_) { observeAndReducer_ = defaultObserveReducerCGM; }
	if (!inplaceDamper_) { inplaceDamper_ = defaultInplaceWeakDamperCGM; }

	// A quick check
	ASSERT( (means.size() == N_) && (covs.size() == N_),
			"weights.size() = " << weights.size() << ", but means.size() = " <<
			means.size() << "and covs.size() = " << covs.size() );

	// Convert from Covariance to Canonical form, this is done upfront
	// as using adjustMass after initialisation is more expensive.
	int fail = 0;
	double detK, g;
	ColVector<double> h;
	Matrix<double> K;

	for (unsigned i = 0; i < N_; i++) {
		K = inv(covs[i], detK, fail);
		if (fail) printf("Could not invert cov[%d] at line number %d in file %s\n", i, __LINE__, __FILE__);

		h = K*means[i];
		g = -0.5*( (K*means[i]).transpose() )*means[i] 
			- log( (pow(2*M_PI, 0.5*vars.size()) / pow(detK, 0.5) ) / weights[i]);

		comps_[i] = uniqptr<Factor> ( new GaussCanonical(vars, K, h, g, false) );
	}

	// Make the sure high level description is sorted.
	vars_ = comps_[0]->getVars();
} // Covariance constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const emdw::RVIds& vars,
		const std::vector<Matrix<double>>& prec,
		const std::vector<ColVector<double>>& info,
		const std::vector<double>& g,
		bool presorted,
		const unsigned maxComponents,
		const double threshold,
		const double unionDistance,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper)
			: vars_(vars.size()),
			comps_(g.size()),
			N_(g.size()),
			maxComp_(maxComponents),
			threshold_(threshold),
			unionDistance_(unionDistance),
			inplaceNormalizer_(inplaceNormalizer),
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
	if (!inplaceNormalizer_) { inplaceNormalizer_ = defaultInplaceNormalizerCGM; }
	if (!normalizer_) { normalizer_ = defaultNormalizerCGM; }
	if (!inplaceAbsorber_) { inplaceAbsorber_ = defaultInplaceAbsorberCGM; }
	if (!absorber_) { absorber_ = defaultAbsorberCGM; }
	if (!inplaceCanceller_) { inplaceCanceller_ = defaultInplaceCancellerCGM; }
	if (!canceller_) { canceller_ = defaultCancellerCGM; }
	if (!marginalizer_) { marginalizer_ = defaultMarginalizerCGM; }
	if (!observeAndReducer_) { observeAndReducer_ = defaultObserveReducerCGM; }
	if (!inplaceDamper_) { inplaceDamper_ = defaultInplaceWeakDamperCGM; }

	// A quick check
	ASSERT( (info.size() == N_) && (prec.size() == N_),
			"g.size() = " << N_ << ", but info.size() = " <<
			info.size() << "and prec.size() = " << prec.size() );

	for (unsigned i = 0; i < N_; i++) {
		comps_[i] = uniqptr<Factor> ( new GaussCanonical(vars, prec[i], info[i], g[i], false ) );
	}

	// Make the sure high level description is sorted.
	vars_ = comps_[0]->getVars();
} // Canonical constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const emdw::RVIds& vars,
		const std::vector<rcptr<Factor>>& components,
		bool presorted,
		const unsigned maxComponents,
		const double threshold,
		const double unionDistance,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) 
			: vars_(vars.size()),
			comps_(components.size()),
			N_(components.size()),
			maxComp_(maxComponents),
			threshold_(threshold),
			unionDistance_(unionDistance),
			inplaceNormalizer_(inplaceNormalizer),
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
	if (!inplaceNormalizer_) { inplaceNormalizer_ = defaultInplaceNormalizerCGM; }
	if (!normalizer_) { normalizer_ = defaultNormalizerCGM; }
	if (!inplaceAbsorber_) { inplaceAbsorber_ = defaultInplaceAbsorberCGM; }
	if (!absorber_) { absorber_ = defaultAbsorberCGM; }
	if (!inplaceCanceller_) { inplaceCanceller_ = defaultInplaceCancellerCGM; }
	if (!canceller_) { canceller_ = defaultCancellerCGM; }
	if (!marginalizer_) { marginalizer_ = defaultMarginalizerCGM; }
	if (!observeAndReducer_) { observeAndReducer_ = defaultObserveReducerCGM; }
	if (!inplaceDamper_) { inplaceDamper_ = defaultInplaceWeakDamperCGM; }
	
	// Make the sure high level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;	
	} else {
		std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
		vars_ = extract<unsigned>(vars, sorted);
	}

	for (unsigned i = 0; i < N_; i++) {
		ASSERT( vars == components[i]->getVars(), vars << " != " << components[i]->getVars()
				<< ". All components must be distributions in " << vars);
		
		comps_[i] = uniqptr<Factor>(components[i]->copy());
	}
} // Component constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const rcptr<Factor>& xFPtr,
		const Matrix<double>& A,
		const emdw::RVIds& newVars,
		const Matrix<double>& Q,
		bool presorted,
		const unsigned maxComponents,
		const double threshold,
		const double unionDistance,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) 
			: maxComp_(maxComponents),
			threshold_(threshold),
			unionDistance_(unionDistance),
			inplaceNormalizer_(inplaceNormalizer),
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
	if (!inplaceNormalizer_) { inplaceNormalizer_ = defaultInplaceNormalizerCGM; }
	if (!normalizer_) { normalizer_ = defaultNormalizerCGM; }
	if (!inplaceAbsorber_) { inplaceAbsorber_ = defaultInplaceAbsorberCGM; }
	if (!absorber_) { absorber_ = defaultAbsorberCGM; }
	if (!inplaceCanceller_) { inplaceCanceller_ = defaultInplaceCancellerCGM; }
	if (!canceller_) { canceller_ = defaultCancellerCGM; }
	if (!marginalizer_) { marginalizer_ = defaultMarginalizerCGM; }
	if (!observeAndReducer_) { observeAndReducer_ = defaultObserveReducerCGM; }
	if (!inplaceDamper_) { inplaceDamper_ = defaultInplaceWeakDamperCGM; }

	// Get the old mixture components.
	rcptr<CanonicalGaussianMixture> cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(xFPtr);
	std::vector<rcptr<Factor>> oldComps = cgm->getComponents();
	
	// Allocate the new mixture.
	N_ = oldComps.size();
	comps_ = std::vector<rcptr<Factor>>(N_);

	// Put each component through the linear transform.
	for (unsigned i = 0; i < N_; i++) {
		comps_[i] = uniqptr<Factor>(new GaussCanonical(oldComps[i].get(), A, newVars, Q, false ) );
	}

	// Make the new variables are sorted in CanonicalGaussianMixture
	vars_ = comps_[0]->getVars();
} // Linear Gaussian constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const rcptr<Factor>& xFPtr,
		const rcptr<V2VTransform>& transform,
		const emdw::RVIds& newVars,
		const Matrix<double>& Q,
		bool presorted,
		const unsigned maxComponents,
		const double threshold,
		const double unionDistance,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& marginalizer,
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) 
			: maxComp_(maxComponents),
			threshold_(threshold),
			unionDistance_(unionDistance),
			inplaceNormalizer_(inplaceNormalizer),
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
	if (!inplaceNormalizer_) { inplaceNormalizer_ = defaultInplaceNormalizerCGM; }
	if (!normalizer_) { normalizer_ = defaultNormalizerCGM; }
	if (!inplaceAbsorber_) { inplaceAbsorber_ = defaultInplaceAbsorberCGM; }
	if (!absorber_) { absorber_ = defaultAbsorberCGM; }
	if (!inplaceCanceller_) { inplaceCanceller_ = defaultInplaceCancellerCGM; }
	if (!canceller_) { canceller_ = defaultCancellerCGM; }
	if (!marginalizer_) { marginalizer_ = defaultMarginalizerCGM; }
	if (!observeAndReducer_) { observeAndReducer_ = defaultObserveReducerCGM; }
	if (!inplaceDamper_) { inplaceDamper_ = defaultInplaceWeakDamperCGM; }

	// Get the old mixture components.
	rcptr<CanonicalGaussianMixture> cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(xFPtr);
	std::vector<rcptr<Factor>> oldComps = cgm->getComponents();
	
	// Allocate the new mixture.
	N_ = oldComps.size();
	comps_ = std::vector<rcptr<Factor>>(N_);

	// Put each component through the transform.
	for (unsigned i = 0; i < N_; i++) {
		comps_[i] = uniqptr<Factor>(new GaussCanonical(oldComps[i].get(), *transform, newVars, Q, false ) );

		// Preserve the mass
		double mass = (std::dynamic_pointer_cast<GaussCanonical>(oldComps[i]))->getMass();
		(std::dynamic_pointer_cast<GaussCanonical>(comps_[i]))->adjustMass(mass);
	}

	// Make the new variables are sorted in CanonicalGaussianMixture
	vars_ = comps_[0]->getVars();	
} // Non-linear Gaussian constructor

CanonicalGaussianMixture::~CanonicalGaussianMixture() {} // Default Destructor

unsigned CanonicalGaussianMixture::configure(unsigned) {
	std::cout << "NIY" << std::endl;
	return true;
} // configure()

unsigned CanonicalGaussianMixture::classSpecificConfigure(
		const emdw::RVIds& vars,
		const std::vector<rcptr<Factor>>& components,
		bool presorted,
		const unsigned maxComponents,
		const double threshold,
		const double unionDistance,
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
	this->~CanonicalGaussianMixture();

	// .. and begin anew!
	new(this) CanonicalGaussianMixture(
			vars,
			components,
			presorted,
			maxComponents,
			threshold,
			unionDistance,
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

inline void CanonicalGaussianMixture::inplaceNormalize(FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this);
	else dynamicInplaceApply(inplaceNormalizer_.get(), this);
} // inplaceNormalize()

inline uniqptr<Factor> CanonicalGaussianMixture::normalize(FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor>(dynamicApply(procPtr, this));
	else return uniqptr<Factor>(dynamicApply(normalizer_.get(), this));
} // normalize()

//------------------Family 2: Absorbtion, Cancellation

inline void CanonicalGaussianMixture::inplaceAbsorb(const Factor* rhsPtr, FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
	else dynamicInplaceApply(inplaceAbsorber_.get(), this, rhsPtr);
} // inplaceAbsorb()

inline uniqptr<Factor> CanonicalGaussianMixture::absorb(const Factor* rhsPtr, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, rhsPtr));
	else return uniqptr<Factor> (dynamicApply(absorber_.get(), this, rhsPtr));
} // absorb()

inline void CanonicalGaussianMixture::inplaceCancel(const Factor* rhsPtr, FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
	else dynamicInplaceApply(inplaceCanceller_.get(), this, rhsPtr);
} // inplaceCancel()

inline uniqptr<Factor> CanonicalGaussianMixture::cancel(const Factor* rhsPtr, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, rhsPtr));
	else return uniqptr<Factor> (dynamicApply(canceller_.get(), this, rhsPtr));
} // cancel()

//------------------Family 4: Marginalization

inline uniqptr<Factor> CanonicalGaussianMixture::marginalize(const emdw::RVIds& variablesToKeep, 
		bool presorted, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, variablesToKeep, presorted));
	else return uniqptr<Factor> (dynamicApply(marginalizer_.get(), this, variablesToKeep, presorted));
} // marginalize()

//------------------Family 4: ObserveAndReduce

inline uniqptr<Factor> CanonicalGaussianMixture::observeAndReduce( const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted, FactorOperator* procPtr) const {
	if (procPtr) return uniqptr<Factor> (dynamicApply(procPtr, this, variables, assignedVals, presorted));
	else return uniqptr<Factor> (dynamicApply(observeAndReducer_.get(), this, variables, assignedVals, presorted));
} // observeAndReduce()


//------------------Family 4: Inplace Weak Damping

// TODO: Complete this!!!
double CanonicalGaussianMixture::inplaceDampen(const Factor* oldMsg, double df, FactorOperator* procPtr) {
	if (procPtr) return dynamicInplaceApply(procPtr, this, oldMsg, df);
	else return dynamicInplaceApply(inplaceDamper_.get(), this, oldMsg, df); 
} // inplaceDampen()

//------------------Other required virtual methods

CanonicalGaussianMixture* CanonicalGaussianMixture::copy(const emdw::RVIds& newVars, bool presorted) const {
	if (newVars.size()) {
		// Copy components onto new scope
		std::vector<rcptr<Factor>> components;
		for (rcptr<Factor> i : comps_) components.push_back( uniqptr<Factor> (i->copy(newVars, presorted) ) );

		return new CanonicalGaussianMixture(
				newVars,
				components,
				true,
				maxComp_,
				threshold_,
				unionDistance_,
				inplaceNormalizer_,
				normalizer_,
				inplaceAbsorber_,
				absorber_,
				inplaceCanceller_,
				canceller_,
				marginalizer_,
				observeAndReducer_,
				inplaceDamper_);
	} 
	return new CanonicalGaussianMixture(*this);
} // copy()

CanonicalGaussianMixture* CanonicalGaussianMixture::vacuousCopy(const emdw::RVIds& selectedVars, bool presorted) const {
	return new CanonicalGaussianMixture(selectedVars);
} // vacuousCopy()

bool CanonicalGaussianMixture::isEqual(const Factor* rhsPtr) const { return true; } // isEqual()

unsigned CanonicalGaussianMixture::noOfVars() const { return vars_.size(); } // noOfVars()

emdw::RVIds CanonicalGaussianMixture::getVars() const { return vars_; } // getVars()

emdw::RVIdType CanonicalGaussianMixture::getVar(unsigned varNo) const { return vars_[varNo]; } // getVar()

//TODO: Complete this!!!
std::istream& CanonicalGaussianMixture::txtRead(std::istream& file) { return file; } // txtRead()

//TODO: Complete this!!
std::ostream& CanonicalGaussianMixture::txtWrite(std::ostream& file) const { 
	std::vector<rcptr<Factor>> components = getComponents();
	
	for (unsigned i = 0; i < N_; i++) {
		file << "\n=========================\n";
		file << "Component " << i << "\n";
		file << *components[i] << "\n\n";
		file << "=========================\n";
	}
	
	return file; 
} // txtWrite()

//------------------ M-Projection

uniqptr<Factor> CanonicalGaussianMixture::momentMatch() const {
	return mProject(comps_);
} // momentMatch()

void CanonicalGaussianMixture::pruneAndMerge() {
	if (N_ > maxComp_) {
		std::vector<rcptr<Factor>> reduced = pruneComponents(comps_, threshold_);
		std::vector<rcptr<Factor>> merged =  mergeComponents(reduced, maxComp_, threshold_, unionDistance_);
		
		emdw::RVIds vars = vars_;
		comps_.clear(); reduced.clear();

		classSpecificConfigure( vars, 
					merged, 
					true,
					maxComp_,
					threshold_,
					unionDistance_,
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
} //pruneAndMerge()

//---------------- Adjust Mass

void CanonicalGaussianMixture::adjustMass(const double mass) {
	for (rcptr<Factor> c : comps_) {
		std::dynamic_pointer_cast<GaussCanonical>(c)->adjustMass(mass);
	}
} // adjustMass()

//---------------- Useful get methods

std::vector<rcptr<Factor>> CanonicalGaussianMixture::getComponents() const { 
	std::vector<rcptr<Factor>> components(comps_.size());

	for (unsigned i = 0; i < comps_.size(); i++) {
		components[i] = uniqptr<Factor>( comps_[i]->copy() ) ;
	} // for
	
	return components; 
} // getComponents()

double CanonicalGaussianMixture::getNumberOfComponents() const { return N_; } // getNumberOfComponents()

std::vector<double> CanonicalGaussianMixture::getWeights() const {
	std::vector<double> weights(comps_.size());
	for (unsigned i = 0; i < comps_.size(); i++) {
		weights[i]  = (std::dynamic_pointer_cast<GaussCanonical>(comps_[i]))->getMass() ;
	} // for
	return weights;
} // getWeights()

std::vector<ColVector<double>> CanonicalGaussianMixture::getMeans() const {
	std::vector<ColVector<double>> means;
	for (rcptr<Factor> c : comps_) means.push_back( (std::dynamic_pointer_cast<GaussCanonical>(c))->getMean() );
	return means;
} // getMeans()

std::vector<Matrix<double>> CanonicalGaussianMixture::getCovs() const {
	std::vector<Matrix<double>> covs;
	for (rcptr<Factor> c : comps_) covs.push_back( (std::dynamic_pointer_cast<GaussCanonical>(c))->getCov() );
	return covs;
} // getCovs()

std::vector<double> CanonicalGaussianMixture::getG() const {
	std::vector<double> g;
	for (rcptr<Factor> c : comps_) g.push_back( (std::dynamic_pointer_cast<GaussCanonical>(c))->getG() );
	return g;
} // getG()

std::vector<ColVector<double>> CanonicalGaussianMixture::getH() const {
	std::vector<ColVector<double>> info;
	for (rcptr<Factor> c : comps_) info.push_back( (std::dynamic_pointer_cast<GaussCanonical>(c))->getH() );
	return info;
} // getH()

std::vector<Matrix<double>> CanonicalGaussianMixture::getK() const {
	std::vector<Matrix<double>> prec;
	for (rcptr<Factor> c : comps_) prec.push_back( (std::dynamic_pointer_cast<GaussCanonical>(c))->getK());
	return prec;
} // getK()

//==================================================FactorOperators======================================

//------------------Family 1: Normalization
//
const std::string& InplaceNormalizeCGM::isA() const {
	static const std::string CLASSNAME("InplaceNormalizeCGM");
	return CLASSNAME;
} // isA()

void InplaceNormalizeCGM::inplaceProcess(CanonicalGaussianMixture* lhsPtr) {
	CanonicalGaussianMixture& lhs(*lhsPtr);
	std::vector<rcptr<Factor>> lhsComp = lhs.getComponents();

	// Get the total mass
	std::vector<double> weights = lhs.getWeights();
	double totalMass = 0.0;
	for (auto& w : weights) totalMass += w;

	// Divide through by the total mass
	// TODO: Sort this out, it is a mess.
	for (rcptr<Factor> c : lhsComp) {
		std::dynamic_pointer_cast<GaussCanonical>(c)->adjustMass(1.0/totalMass);
	}

	// Reconfigure
	lhs.classSpecificConfigure(lhs.getVars(), 
			lhsComp, 
			true, 
			lhs.maxComp_,
			lhs.threshold_,
			lhs.unionDistance_,
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

const std::string& NormalizeCGM::isA() const {
	static const std::string CLASSNAME("NormalizeCGM");
	return CLASSNAME;
} // isA()

Factor* NormalizeCGM::process(const CanonicalGaussianMixture* lhsPtr) {
	CanonicalGaussianMixture* fPtr = new CanonicalGaussianMixture(*lhsPtr);
	InplaceNormalizeCGM ipNorm;
	
	try { 
		ipNorm.inplaceProcess(fPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()


//------------------Family 2: Absorption, Cancellation

const std::string& InplaceAbsorbCGM::isA() const {
	static const std::string CLASSNAME("InplaceAbsorbCGM");
	return CLASSNAME;
} // isA()

void InplaceAbsorbCGM::inplaceProcess(CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {
	// Try cast the pointer to CanonicalGaussianMixture 
	CanonicalGaussianMixture& lhs(*lhsPtr);
	const CanonicalGaussianMixture* rhsCGMPtr = dynamic_cast<const CanonicalGaussianMixture*>(rhsFPtr);
	
	// New components
	std::vector<rcptr<Factor>> product;
	std::vector<rcptr<Factor>> lhsComps = lhs.getComponents();

	// If it isn't a CanonicalGaussianMixture then GaussCanonical should do all the validation.
	if (rhsCGMPtr){
		const CanonicalGaussianMixture& rhs(*rhsCGMPtr);
		std::vector<rcptr<Factor>> rhsComps = rhs.getComponents();

		for (rcptr<Factor> i : lhsComps) {
			for (rcptr<Factor> j : rhsComps) {
				product.push_back(i->absorb(j));
			}
		}
	} else { 
		rcptr<Factor> rhs = uniqptr<Factor>(rhsFPtr->copy());
		for (rcptr<Factor> c : lhsComps) product.push_back(c->absorb(rhs));
	}
	
	emdw::RVIds vars = product[0]->getVars();

	// Reconfigure the class
	lhs.classSpecificConfigure(vars, 
				product, 
				true,
				lhs.maxComp_,
				lhs.threshold_,
				lhs.unionDistance_,
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

const std::string& AbsorbCGM::isA() const {
	static const std::string CLASSNAME("AbsorbCGM");
	return CLASSNAME;
} // isA()

Factor* AbsorbCGM::process(const CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {
	CanonicalGaussianMixture* fPtr = new CanonicalGaussianMixture(*lhsPtr);
	InplaceAbsorbCGM ipAbsorb;
	
	try { 
		ipAbsorb.inplaceProcess(fPtr, rhsFPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()

const std::string& InplaceCancelCGM::isA() const {
	static const std::string CLASSNAME("InplaceCancelCGM");
	return CLASSNAME;
} // isA()

void InplaceCancelCGM::inplaceProcess(CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {
	// Try cast the pointer to CanonicalGaussianMixture 
	CanonicalGaussianMixture& lhs(*lhsPtr);
	const CanonicalGaussianMixture* rhsCGMPtr = dynamic_cast<const CanonicalGaussianMixture*>(rhsFPtr);
	const CanonicalGaussianMixture& rhs(*rhsCGMPtr); // Not very safe
	
	// New components
	std::vector<rcptr<Factor>> lhsComps = lhs.getComponents();
	rcptr<Factor> single;
	
	// If it isn't a CanonicalGaussianMixture then GaussCanonical should do all the validation.
	if (rhsCGMPtr) single = rhs.momentMatch();
	else single = uniqptr<Factor>(rhsFPtr->copy());

	// Divide through by a single GaussCanonical
	std::vector<rcptr<Factor>> quotient; quotient.clear();
	for (rcptr<Factor> i : lhsComps) { 
		quotient.push_back( i->cancel(single)  ); 

		//double mass = std::dynamic_pointer_cast<GaussCanonical>(quotient.back())->getMass();
		//std::cout << "InplaceCancelCGM, mass: " << mass << std::endl;
	} // for

	// Reconfigure the class
	lhs.classSpecificConfigure(quotient[0]->getVars(), quotient, true,
				lhs.maxComp_,
				lhs.threshold_,
				lhs.unionDistance_,
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

const std::string& CancelCGM::isA() const {
	static const std::string CLASSNAME("CancelCGM");
	return CLASSNAME;
} // isA()


Factor* CancelCGM::process(const CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {
	CanonicalGaussianMixture* fPtr = new CanonicalGaussianMixture(*lhsPtr);
	InplaceCancelCGM ipCancel;
	
	try { 
		ipCancel.inplaceProcess(fPtr, rhsFPtr); 
	} catch (const char* s) {
		std::cout << __FILE__ << __LINE__ << " call to 'inplaceProcess' failed" << std::endl;
		throw s;
	}

	return fPtr;
} // process()


//------------------Family 3: Marginalization

const std::string& MarginalizeCGM::isA() const {
	static const std::string CLASSNAME("MarginalizeCGM");
	return CLASSNAME;
} // isA()

Factor* MarginalizeCGM::process(const CanonicalGaussianMixture* lhsPtr, const emdw::RVIds& variablesToKeep,
		bool presorted) {
	const CanonicalGaussianMixture& lhs(*lhsPtr);
	std::vector<rcptr<Factor>> lhsComps = lhs.getComponents();
	unsigned M = lhsComps.size();	
	std::vector<rcptr<Factor>> result(M);

	// If everything is marginalized out.
	if (!variablesToKeep.size()) {
		return new CanonicalGaussianMixture(
				variablesToKeep,
				true);
	}

	// Let GaussCanonical sort it all out for us.
	for (unsigned i = 0; i < M; i++) {
		result[i] = (lhsComps[i])->marginalize(variablesToKeep, presorted);
		result[i]->inplaceNormalize();

		double weight = (std::dynamic_pointer_cast<GaussCanonical>(lhsComps[i]))->getMass();
		(std::dynamic_pointer_cast<GaussCanonical>(result[i])->adjustMass(weight));
	}

	return new CanonicalGaussianMixture(result[0]->getVars(), 
				result, 
				true,
				lhs.maxComp_,
				lhs.threshold_,
				lhs.unionDistance_,
				lhs.inplaceNormalizer_,
				lhs.normalizer_,
				lhs.inplaceAbsorber_,
				lhs.absorber_,
				lhs.inplaceCanceller_,
				lhs.canceller_,
				lhs.marginalizer_,
				lhs.observeAndReducer_,
				lhs.inplaceDamper_ );
} // process()


//------------------Family 4: ObserveAndReduce

const std::string& ObserveAndReduceCGM::isA() const {
	static const std::string CLASSNAME("ObserveAndReduceCGM");
	return CLASSNAME;
} // isA()

Factor* ObserveAndReduceCGM::process(const CanonicalGaussianMixture* lhsPtr, const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted) {
	const CanonicalGaussianMixture& lhs(*lhsPtr);
	std::vector<rcptr<Factor>> lhsComps = lhs.getComponents();
	std::vector<rcptr<Factor>> result;

	// If nothing was observed.
	if(!variables.size()) return lhs.copy(); 

	// Let GaussCanonical sort it all out for us.
	for (rcptr<Factor> c : lhsComps) result.push_back(c->observeAndReduce(variables, assignedVals, presorted));

	//for (rcptr<Factor> c : result) std::cout << "CGM.observeAndReduce(): " << std::dynamic_pointer_cast<GaussCanonical>(c)->getMass() << std::endl;

	return new CanonicalGaussianMixture(result[0]->getVars(), 
			        result, 
				true,
				lhs.maxComp_,
				lhs.threshold_,
				lhs.unionDistance_,
				lhs.inplaceNormalizer_,
				lhs.normalizer_,
				lhs.inplaceAbsorber_,
				lhs.absorber_,
				lhs.inplaceCanceller_,
				lhs.canceller_,
				lhs.marginalizer_,
				lhs.observeAndReducer_,
				lhs.inplaceDamper_ );
} // process()


//------------------Family 5: Damping

const std::string& InplaceWeakDampingCGM::isA() const {
	static const std::string CLASSNAME("InplaceWeakDampingCGM");
	return CLASSNAME;
} // isA()

// TODO: Complete this!!!
double InplaceWeakDampingCGM::inplaceProcess(const CanonicalGaussianMixture* lhsPtr, const Factor* rhsPtr, double df) {
	return 0.0;
} // inplaceProcess()

//------------------ M-Projections

uniqptr<Factor> mProject(const std::vector<rcptr<Factor>>& components) {
	// Old means and covariances
	unsigned M = components.size();
	ASSERT( M != 0, "There must be at least one mixand" );

	// The GM's components
	std::vector<double> w(M);
	std::vector<ColVector<double>> mu(M);
	std::vector<Matrix<double>> S(M);
	double totalMass = 0.0;

	// Scope and dimension
	emdw::RVIds vars = (components.back())->getVars();
	unsigned dimension = vars.size();

	// First and second central moments
	ColVector<double> mean(dimension); mean *= 0.0;
	Matrix<double> cov = gLinear::zeros<double>(dimension, dimension); cov *= 0.0;

	// Get the non-vacuous Gaussians' weights, means and covariances - Not very efficient, but whatever.
	for (unsigned i = 0; i < M; i++) {
		if (components[i]->noOfVars() != 0) { // Horrible hack to check if vacuous
			rcptr<Factor> factor = uniqptr<Factor>( components[i]->copy()  );
			rcptr<GaussCanonical> gc = std::dynamic_pointer_cast<GaussCanonical>(factor);

			// Get the mean and covariance
			mu[i] = gc->getMean();
			S[i] = gc->getCov();

			// Determine the weight
			w[i] = gc->getMass();
			totalMass += w[i];
		} else {
			w[i] = 0.0; 
			mu[i] = ColVector<double>(dimension); mu[i] *= 0;
			S[i] = gLinear::zeros<double>(dimension, dimension); S[i] *= 0;
		} // if
	} // for

	// Determine the first two moments of the mixture
	for (unsigned i = 0; i < M; i++) {
		double weight = w[i]/totalMass;
		mean += (weight)*(mu[i]);
		cov += (weight)*( S[i] + (mu[i])*(mu[i].transpose())  );
	} // for
	cov -= (mean)*(mean.transpose());

	return uniqptr<Factor>(new GaussCanonical(vars, mean, cov) );
} // mProject()

//------------------ Pruning and Merging

std::vector<rcptr<Factor>> pruneComponents(const std::vector<rcptr<Factor>>& components, const double threshold)  {
	std::vector<rcptr<Factor>> reduced; reduced.clear();

	for (rcptr<Factor> c : components) {
		double mass = std::dynamic_pointer_cast<GaussCanonical>(c)->getMass();
		if (mass > threshold) reduced.push_back( uniqptr<Factor>(c->copy()) );
	} // for
	return reduced;
} // pruneComponents()

std::vector<rcptr<Factor>> mergeComponents(const std::vector<rcptr<Factor>>& components, const unsigned maxComp,
		const double threshold, const double unionDistance) {
	ASSERT( components.size() != 0, "There must be at least one mixand." );

	/*
	std::cout << "Before merge\n" << std::endl;
	for (unsigned i = 0; i < components.size(); i++) {
		rcptr<GaussCanonical> gc = std::dynamic_pointer_cast<GaussCanonical>(components[i]);
		double mass = gc->getMass();
		ColVector<double> mean = gc->getMean();
		Matrix<double> cov = gc->getCov();

		std::cout << "Component " << i << std::endl;
		std::cout << "Mass: " << mass << std::endl;
		std::cout << "Mean: " << mean << std::endl;
		std::cout << "Cov: " << cov << std::endl;
	}
	*/

	//std::cout << "Number of components: " << components.size() << std::endl;

	// Local variables
	emdw::RVIds vars = (components.back())->getVars();
	std::vector<rcptr<Factor>> merged; merged.clear();
	std::vector<rcptr<Factor>> comps; comps.clear();
	std::vector<double> weights; weights.clear();

	// Get the weights and copy the factors
	for (rcptr<Factor> c : components) {
		comps.push_back( uniqptr<Factor>( c->copy() ) );
		weights.push_back( std::dynamic_pointer_cast<GaussCanonical>(c)->getMass()  );
	} // for	

	// Sort the components according to weight
	std::vector<size_t> sortedIndices = sortIndices( weights, std::less<double>() );
	std::vector<double> w = extract<double>( weights, sortedIndices );
	std::vector<rcptr<Factor>> newComps = extract<rcptr<Factor>>(comps, sortedIndices);

	// Merge closely spaced components
	while ( !newComps.empty() ) {
		// Dominant components mean
		ColVector<double> mu_0 = std::dynamic_pointer_cast<GaussCanonical>(newComps[0])->getMean();
		unsigned L = w.size();

		// Local variables
		std::map<unsigned, bool> indices; indices.clear();
		ColVector<double> mu( vars.size() ); mu *= 0;
		Matrix<double> S = gLinear::zeros<double>( vars.size(), vars.size() );
		double g = 0;

		// Create a merged super Gaussian
		for (unsigned i = 0; i < L; i++) {
			rcptr<GaussCanonical> gc = std::dynamic_pointer_cast<GaussCanonical>( newComps[i]  );
			if (gc->mahalanobis(mu_0) <= unionDistance) {
				//std::cout << "Merged: " << i << std::endl;
				ColVector<double> mean = gc->getMean();
				
				g += w[i];
				mu += w[i]*(mean);
				S += w[i]*( gc->getCov() + (mean - mu_0)*( (mean - mu_0).transpose() ) );
				
				indices[i] = true;
			} // if
		} // for

		// Remove merged indices from the list -- Terrible method, but safe-ish.
		std::vector<rcptr<Factor>> cTemp; cTemp.clear();
		std::vector<double> wTemp; wTemp.clear();
		for (unsigned i = 0; i < L; i++) {
			if (!indices[i]) {
				cTemp.push_back(newComps[i]);
				wTemp.push_back(w[i]);
			} // if
		} // for

		// Reassign temporary variables to local variables
		newComps.clear(); w.clear();
		for (unsigned i = 0 ; i < cTemp.size(); i++) {
			newComps.push_back(cTemp[i]);
			w.push_back(wTemp[i]);
		} // for

		// Create a new Gaussian
		rcptr<Factor> scaled = uniqptr<Factor> (new GaussCanonical( vars, mu/g, S/g, true ));
		std::dynamic_pointer_cast<GaussCanonical>(scaled)->adjustMass(g);
		merged.push_back(scaled);
	} // while	

	//std::cout << "Number of components: " << merged.size() << std::endl;

	// If there are still too many components
	if (merged.size() > maxComp) {
		// Local variables
		std::vector<rcptr<Factor>> sigComps; sigComps.clear();
		std::vector<double> sigWeights; sigWeights.clear();
		
		// Extract weights and copy
		for ( rcptr<Factor> c : merged ) {
			sigComps.push_back( uniqptr<Factor> (c->copy() ) );
			sigWeights.push_back( std::dynamic_pointer_cast<GaussCanonical>(c)->getMass() );
		} // for
		
		// Sort according to weight
		std::vector<size_t> sorted = sortIndices( sigWeights, std::less<double>() );
		std::vector<rcptr<Factor>> ordered = extract<rcptr<Factor>>(sigComps, sorted);
		
		// Select only the N largest components	
		merged.clear();
		for (unsigned i = 0; i < maxComp; i++) merged.push_back(ordered[i]);
	} // if

	/*
	std::cout << "After merge\n" << std::endl;
	for (unsigned i = 0; i < merged.size(); i++) {
		rcptr<GaussCanonical> gc = std::dynamic_pointer_cast<GaussCanonical>(merged[i]);
		double mass = gc->getMass();
		ColVector<double> mean = gc->getMean();
		Matrix<double> cov = gc->getCov();

		std::cout << "Component " << i << std::endl;
		std::cout << "Mass: " << mass << std::endl;
		std::cout << "Mean: " << mean << std::endl;
		std::cout << "Cov: " << cov << std::endl;
	}
	*/

	return merged;
} // mergeComponents()

