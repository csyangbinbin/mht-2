/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for the Canonical Gaussian Mixture implementation.
 *************************************************************************/
#include <vector>
#include <iostream>
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
			unionDistance_(unionDistance)
		{
	
	// Default operator intialisation
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerCGM;
	if (!normalizer) normalizer_ = defaultNormalizerCGM;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorberCGM;
	if (!absorber) absorber_ = defaultAbsorberCGM;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCancellerCGM;
	if (!canceller) canceller_ = defaultCancellerCGM;
	if (!marginalizer) marginalizer_ = defaultMarginalizerCGM;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducerCGM;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamperCGM;

	// Ensure the higher level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;
	} else {
		std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
		vars_ = extract<unsigned>(vars, sorted);
	}

	// Create a mixture with a single vacuous component.
	comps_.push_back( uniqptr<GaussCanonical> ( new GaussCanonical(vars_, true ) ) );
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
			unionDistance_(unionDistance)
		{

	// Default operator initialisation
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerCGM;
	if (!normalizer) normalizer_ = defaultNormalizerCGM;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorberCGM;
	if (!absorber) absorber_ = defaultAbsorberCGM;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCancellerCGM;
	if (!canceller) canceller_ = defaultCancellerCGM;
	if (!marginalizer) marginalizer_ = defaultMarginalizerCGM;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducerCGM;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamperCGM;

	// A quick check
	ASSERT( (means.size() == N_) && (covs.size() == N_),
			"weights.size() = " << weights.size() << ", but means.size() = " <<
			means.size() << "and covs.size() = " << covs.size() );

	// Convert from Covariance to Canonical form, this is done upfront
	// as using adjustMass after initialisation is more expensive.
	int fail;
	double detcov, g;
	ColVector<double> h;
	Matrix<double> K;

	for (unsigned i = 0; i < N_; i++) {
		K = inv(covs[i], detcov, fail);
		if (fail) printf("Could not invert cov[%d] at line number %d in file %s\n", i, __LINE__, __FILE__);

		h = K*means[i];
		g = -0.5*( (K*means[i]).transpose() )*means[i] 
			- log( ( pow(2*M_PI, (1.0*vars.size())/2) * pow(detcov, 0.5) ) / weights[i]);

		comps_[i] = uniqptr<GaussCanonical> ( new GaussCanonical(vars, K, h, g, false) );
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
			unionDistance_(unionDistance)
		{

	// Default initialisation
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerCGM;
	if (!normalizer) normalizer_ = defaultNormalizerCGM;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorberCGM;
	if (!absorber) absorber_ = defaultAbsorberCGM;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCancellerCGM;
	if (!canceller) canceller_ = defaultCancellerCGM;
	if (!marginalizer) marginalizer_ = defaultMarginalizerCGM;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducerCGM;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamperCGM;

	// A quick check
	ASSERT( (info.size() == N_) && (prec.size() == N_),
			"g.size() = " << N_ << ", but info.size() = " <<
			info.size() << "and prec.size() = " << prec.size() );

	for (unsigned i = 0; i < N_; i++) {
		comps_[i] = uniqptr<GaussCanonical> ( new GaussCanonical(vars, prec[i], info[i], g[i], false ) );
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
			unionDistance_(unionDistance)
		{
	
	// Default initialisation	
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerCGM;
	if (!normalizer) normalizer_ = defaultNormalizerCGM;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorberCGM;
	if (!absorber) absorber_ = defaultAbsorberCGM;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCancellerCGM;
	if (!canceller) canceller_ = defaultCancellerCGM;
	if (!marginalizer) marginalizer_ = defaultMarginalizerCGM;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducerCGM;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamperCGM;
	
	// Make the sure high level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;	
	} else {
		std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
		vars_ = extract<unsigned>(vars, sorted);
	}

	rcptr<GaussCanonical> gc;
	for (unsigned i = 0; i < N_; i++) {
		ASSERT( vars == components[i]->getVars(), vars << " != " << components[i]->getVars()
				<< ". All components must be distributions in " << vars);
		
		gc = std::dynamic_pointer_cast<GaussCanonical>(components[i]);
		comps_[i] = uniqptr<GaussCanonical>(gc->copy());
	}

	// Prune and merge if necessary - or something
	if (N_ > maxComp_) {
		pruneComponents();
		mergeComponents();
	}
} // Component constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const rcptr<Factor> xFPtr,
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
		) {

	// Default initialisation.
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerCGM;
	if (!normalizer) normalizer_ = defaultNormalizerCGM;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorberCGM;
	if (!absorber) absorber_ = defaultAbsorberCGM;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCancellerCGM;
	if (!canceller) canceller_ = defaultCancellerCGM;
	if (!marginalizer) marginalizer_ = defaultMarginalizerCGM;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducerCGM;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamperCGM;

	// Get the old mixture components.
	rcptr<CanonicalGaussianMixture> cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(xFPtr);
	std::vector<rcptr<Factor>> oldComps = cgm->getComponents();
	
	// Allocate the new mixture.
	N_ = oldComps.size();
	comps_ = std::vector<rcptr<Factor>>(N_);

	// Put each component through the linear transform.
	for (unsigned i = 0; i < N_; i++) {
		comps_[i] = uniqptr<GaussCanonical>(new GaussCanonical(oldComps[i]->copy(), A, newVars, Q, false ) );
	}

	// Make the new variables are sorted in CanonicalGaussianMixture
	vars_ = comps_[0]->getVars();
} // Linear Gaussian constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const rcptr<Factor> xFPtr,
		const rcptr<V2VTransform> transform,
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
		) {

	// Default initialisation.
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizerCGM;
	if (!normalizer) normalizer_ = defaultNormalizerCGM;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorberCGM;
	if (!absorber) absorber_ = defaultAbsorberCGM;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCancellerCGM;
	if (!canceller) canceller_ = defaultCancellerCGM;
	if (!marginalizer) marginalizer_ = defaultMarginalizerCGM;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducerCGM;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamperCGM;

	// Get the old mixture components.
	rcptr<CanonicalGaussianMixture> cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(xFPtr);
	std::vector<rcptr<Factor>> oldComps = cgm->getComponents();
	
	// Allocate the new mixture.
	N_ = oldComps.size();
	comps_ = std::vector<rcptr<Factor>>(N_);

	// Put each component through the transform.
	for (unsigned i = 0; i < N_; i++) {
		comps_[i] = uniqptr<GaussCanonical>(new GaussCanonical(oldComps[i]->copy(), *transform, newVars, Q, false ) );
	}

	// Make the new variables are sorted in CanonicalGaussianMixture
	vars_ = comps_[0]->getVars();
} // Non-linear Gaussian constructor

CanonicalGaussianMixture::~CanonicalGaussianMixture() {} // Default Destructor


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
	return new CanonicalGaussianMixture();
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
		file << "=========================\n";
		file << "Component " << i << "\n";
		file << *components[i] << "\n\n";
		file << "=========================\n";
	}
	
	return file; 
} // txtWrite()

//------------------ M-Projections

uniqptr<Factor> CanonicalGaussianMixture::momentMatch() const {
	// Old means and covariances
	std::vector<double> w;
	std::vector<ColVector<double>> mu;
	std::vector<Matrix<double>> S;
	double totalMass = 0;

	// New mean and covariance
	unsigned dimension = vars_.size();
	ColVector<double> mean(dimension); mean *= 0;
	Matrix<double> cov = gLinear::zeros<double>(dimension, dimension); 

	// Get the non-vacuous Gaussians' weights, means and covariances
	rcptr<GaussCanonical> gc;
	for (rcptr<Factor> c : comps_) {
		if (c->noOfVars() != 0) { // Horrible hack to check if vacuous
			gc = std::dynamic_pointer_cast<GaussCanonical>(c);
			w.push_back( gc->getMass() );
			mu.push_back( gc->getMean() );
			S.push_back( gc->getCov() );
			totalMass += w.back();
		}
	}

	// Determine the first two moments of the mixture
	double weight;
	for (unsigned i = 0; i < w.size(); i++) {
		weight = w[i]/totalMass;
		mean += (weight)*(mu[i]);
		cov += (weight)*( S[i] + (mu[i])*(mu[i].transpose())  );
	}
	cov -= (mean)*(mean.transpose());

	return uniqptr<Factor>(new GaussCanonical(vars_, mean, cov) );
}

//------------------ Pruning and Merging

void CanonicalGaussianMixture::pruneComponents()  {
	std::vector<double> weights = this->getWeights();
	std::vector<rcptr<Factor>> newComps;

	// Throw out all insignificant components
	for (unsigned i = 0; i < N_; i++) {
		if (weights[i] > threshold_) newComps.push_back(comps_[i]); 
	}

	// Re-assign
	comps_ = newComps;
	N_ = comps_.size();
} // pruneComponents()

void CanonicalGaussianMixture::mergeComponents() {
	std::vector<rcptr<Factor>> merged;
	std::vector<rcptr<GaussCanonical>> oldComps, cTemp;
	
	// Weights of old components
	std::vector<double> w, wTemp;
	std::set<unsigned> indices;
	unsigned L, maxIndex;

	// Parameters of new components
	unsigned dim = vars_.size();
	ColVector<double> mu(dim), mu_0(dim), mean(dim);
	Matrix<double> S = gLinear::zeros<double>(dim, dim);
	double g;

	// Get the weights and cast the Factors down to CG
	for (rcptr<Factor> c : comps_) {
		oldComps.push_back(std::dynamic_pointer_cast<GaussCanonical>(c));
		w.push_back( (oldComps.back())->getMass() );
	}

	// Calculate the merged components
	while (!oldComps.empty()) {
		// Determine maximum remaining weight.
		maxIndex = std::max_element(w.begin(), w.end()) - w.begin(); // Hack to get max element's index.
		mu_0 = oldComps[maxIndex]->getMean(); // Dominant component's mean.
		L = w.size();
		
		// Determined merged components.
		g = 0; mu *= 0; mean *= 0; S *= 0;
		for (unsigned i = 0; i < L; i++) {
			if (oldComps[i]->mahanalobis(mu_0) <= unionDistance_) {
				mean = oldComps[i]->getMean();
				
				g += w[i];
				mu += w[i]*(mean);
				S +=  w[i]*(oldComps[i]->getCov() + (mean - mu_0)*((mean - mu_0).transpose()));

				indices.insert(i);
			}
		}

		// Absolutely terrible removal technique, but neat iterator solutions hate rcptr<Factor>.
		for (unsigned i = 0; i < L; i++) {
			if (!indices.count(i)) {
				cTemp.push_back(oldComps[i]);
				wTemp.push_back(w[i]);
			}		
		}
		
		// Re-assign and clear temporary variables
		oldComps = cTemp; cTemp.clear(); 
		w = wTemp; wTemp.clear();
		indices.clear();

		// Create a new Gaussian
		merged.push_back( uniqptr<GaussCanonical> (new GaussCanonical(vars_, mu/g, S/g, true)) );
	}

	// Re-assign
	comps_ = merged;
	N_ = comps_.size();

	// If there are still too many components, throw out the smallest ones
	if (N_ > maxComp_) {
		comps_.clear(); oldComps.clear(); w.clear();

		// Get the new weights
		for (rcptr<Factor> c : merged) w.push_back( (std::dynamic_pointer_cast<GaussCanonical>(c))->getMass() );

		// Sort
		std::vector<size_t> sorted = sortIndices( w, std::less<double>() );
		merged = extract<rcptr<Factor>>(merged, sorted);
		
		// Re-assign
		std::copy(merged.begin(), merged.begin() + maxComp_, std::back_inserter(comps_));
		N_ = maxComp_;
	}
} // mergeComponents()


//---------------- Useful get methods

std::vector<rcptr<Factor>> CanonicalGaussianMixture::getComponents() const { 
	std::vector<rcptr<Factor>> components;

	for (rcptr<Factor> c : comps_) components.push_back(uniqptr<Factor>( c->copy() ) );
	
	return components; 
} // getComponents()

double CanonicalGaussianMixture::getNumberOfComponents() const { return N_; } // getNumberOfComponents()

std::vector<double> CanonicalGaussianMixture::getWeights() const {
	std::vector<double> weights;
	for (rcptr<Factor> c : comps_) weights.push_back( (std::dynamic_pointer_cast<GaussCanonical>(c))->getMass() );
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
	double totalMass = 0;
	for (auto& w : weights) totalMass += w;

	// Divide through by the total mass
	// TODO: Sort this out, it is a mess.
	rcptr<GaussCanonical> gc;
	for (rcptr<Factor> c : lhsComp) {
		gc = std::dynamic_pointer_cast<GaussCanonical>(c);
		gc->classSpecificConfigure(gc->getVars(), gc->getK(), gc->getH(),
				gc->getG() - log(totalMass), true);
	}

	// Reconfigure
	lhs.classSpecificConfigure(lhs.getVars(), lhsComp, true,
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

	// Reconfigure the class
	lhs.classSpecificConfigure(product[0]->getVars(), product, true,
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
	std::vector<rcptr<Factor>> quotient;
	std::vector<rcptr<Factor>> lhsComps = lhs.getComponents();
	rcptr<Factor> single;
	
	// If it isn't a CanonicalGaussianMixture then GaussCanonical should do all the validation.
	if (rhsCGMPtr) single = rhs.momentMatch();
	else single = uniqptr<Factor>(rhsFPtr->copy());

	// Divide through by a single GaussCanonical
	for (rcptr<Factor> i : lhsComps) quotient.push_back( i->cancel(single)  );

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
	InplaceAbsorbCGM ipCancel;
	
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
	std::vector<rcptr<Factor>> result;

	// If everything is marginalized out.
	if (!variablesToKeep.size()) {
		return new CanonicalGaussianMixture(
				variablesToKeep,
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
				lhs.inplaceDamper_
				);
	}

	// Let GaussCanonical sort it all out for us.
	for (rcptr<Factor> c : lhsComps) result.push_back(c->marginalize(variablesToKeep, presorted));

	return new CanonicalGaussianMixture(result[0]->getVars(), result, true,
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

	return new CanonicalGaussianMixture(result[0]->getVars(), result, true,
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
