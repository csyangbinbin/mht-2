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
rcptr<FactorOperator> defaultInplaceNormalizer = uniqptr<FactorOperator>(new InplaceNormalizeCGM());
rcptr<FactorOperator> defaultNormalizer = uniqptr<FactorOperator>(new NormalizeCGM());
rcptr<FactorOperator> defaultInplaceAbsorber = uniqptr<FactorOperator>(new InplaceAbsorbCGM());
rcptr<FactorOperator> defaultAbsorber = uniqptr<FactorOperator>(new AbsorbCGM());
rcptr<FactorOperator> defaultInplaceCanceller = uniqptr<FactorOperator>(new InplaceCancelCGM());
rcptr<FactorOperator> defaultCanceller = uniqptr<FactorOperator>(new CancelCGM());
rcptr<FactorOperator> defaultObserveReducer = uniqptr<FactorOperator>(new ObserveAndReduceCGM());
rcptr<FactorOperator> defaultInplaceWeakDamper = uniqptr<FactorOperator>(new InplaceWeakDampingCGM());

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
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper)
			: vars_(vars.size()),
			maxComp_(maxComponents),
			threshold_(threshold),
			unionDistance_(unionDistance)
		{
	
	// Default operator intialisation
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizer;
	if (!normalizer) normalizer_ = defaultNormalizer;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorber;
	if (!absorber) absorber_ = defaultAbsorber;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCanceller;
	if (!canceller) canceller_ = defaultCanceller;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducer;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamper;

	// Ensure the higher level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;
	} else {
		std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
		vars_ = extract<unsigned>(vars, sorted);
	}

	// Create a mixture with a single vacuous component.
	comps_.push_back( uniqptr<GaussCanonical> ( new GaussCanonical(vars_, true,
				inplaceNormalizer_,
				normalizer_,
				inplaceAbsorber_,
				absorber_,
				inplaceCanceller_,
				canceller_,
				observeAndReducer_,
				inplaceDamper_ ) ) );
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
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizer;
	if (!normalizer) normalizer_ = defaultNormalizer;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorber;
	if (!absorber) absorber_ = defaultAbsorber;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCanceller;
	if (!canceller) canceller_ = defaultCanceller;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducer;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamper;

	// A quick check
	ASSERT( (means.size() == N_) && (covs.size() == N_),
			"weights.size() = " << weights.size() << ", but means.size() = " <<
			means.size() << "and covs.size() = " << covs.size() );

	// Convert from Covariance to Canonical form, the alternative is to 
	// create a covariance GaussCanonical and use adjustMass and dynamic pointer casts.
	// NOTE: adjustMass currently has an integer arithmetic error.
	int fail;
	double detcov;
	double g;
	ColVector<double> h;
	Matrix<double> K;

	for (unsigned i = 0; i < N_; i++) {
		K = inv(covs[i], detcov, fail);
		if (fail) printf("Could not invert cov[%d] at line number %d in file %s\n", i, __LINE__, __FILE__);

		h = K*means[i];
		g = -0.5*( (K*means[i]).transpose() )*means[i] 
			- log( ( pow(2*M_PI, (1.0*vars.size())/2) * pow(detcov, 0.5) ) / weights[i]);

		comps_[i] = uniqptr<GaussCanonical> ( new GaussCanonical(vars, K, h, g,
					false,
					inplaceNormalizer_,
					normalizer_,
					inplaceAbsorber_,
					absorber_,
					inplaceCanceller_,
					canceller_,
					observeAndReducer_,
					inplaceDamper_ ) );
	}

	// Make the sure high level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;
		return;
	} 
	
	std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
	vars_ = extract<unsigned>(vars, sorted);
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
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizer;
	if (!normalizer) normalizer_ = defaultNormalizer;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorber;
	if (!absorber) absorber_ = defaultAbsorber;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCanceller;
	if (!canceller) canceller_ = defaultCanceller;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducer;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamper;

	// A quick check
	ASSERT( (info.size() == N_) && (prec.size() == N_),
			"g.size() = " << N_ << ", but info.size() = " <<
			info.size() << "and prec.size() = " << prec.size() );

	for (unsigned i = 0; i < N_; i++) {
		comps_[i] = uniqptr<GaussCanonical> ( new GaussCanonical(vars, prec[i], info[i], g[i],
					false,
					inplaceNormalizer_,
					normalizer_,
					inplaceAbsorber_,
					absorber_,
					inplaceCanceller_,
					canceller_,
					observeAndReducer_,
					inplaceDamper_ ) );
	}

	// Make the sure high level description is sorted.
	if (presorted || !vars.size()) {
		vars_ = vars;
		return;
	} 
	
	std::vector<size_t> sorted = sortIndices(vars, std::less<unsigned>() );
	vars_ = extract<unsigned>(vars, sorted);
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
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizer;
	if (!normalizer) normalizer_ = defaultNormalizer;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorber;
	if (!absorber) absorber_ = defaultAbsorber;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCanceller;
	if (!canceller) canceller_ = defaultCanceller;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducer;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamper;
	
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

	// Prune and merge if necessary
	pruneComponents();
	mergeComponents();
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
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {

	// Default initialisation
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizer;
	if (!normalizer) normalizer_ = defaultNormalizer;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorber;
	if (!absorber) absorber_ = defaultAbsorber;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCanceller;
	if (!canceller) canceller_ = defaultCanceller;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducer;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamper;

	// Get the old components
	rcptr<CanonicalGaussianMixture> cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(xFPtr);
	std::vector<rcptr<Factor>> oldComps = cgm->getComponents();
	
	// Put each component through the linear transform
	unsigned N = oldComps.size();
	comps_ = std::vector<rcptr<Factor>>(N);

	for (unsigned i = 0; i < N; i++) {
		comps_[i] = uniqptr<GaussCanonical>(new GaussCanonical(oldComps[i]->copy(), A, newVars,Q, false,
					inplaceNormalizer_,
					normalizer_,
					inplaceAbsorber_,
					absorber_,
					inplaceCanceller_,
					canceller_,
					observeAndReducer_,
					inplaceDamper_ ) );
	}

} // Linear Gaussian constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const Factor* xFPtr,
		const V2VTransform& transform,
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
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizer;
	if (!normalizer) normalizer_ = defaultNormalizer;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorber;
	if (!absorber) absorber_ = defaultAbsorber;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCanceller;
	if (!canceller) canceller_ = defaultCanceller;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducer;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamper;
} // Non-linear Gaussian constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const Factor* x1FPtr,
		const Factor* x2FPtr,
		const V2VTransform& transform,
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
		const rcptr<FactorOperator>& observerAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {
	if (!inplaceNormalizer) inplaceNormalizer_ = defaultInplaceNormalizer;
	if (!normalizer) normalizer_ = defaultNormalizer;
	if (!inplaceAbsorber) inplaceAbsorber_ = defaultInplaceAbsorber;
	if (!absorber) absorber_ = defaultAbsorber;
	if (!inplaceCanceller) inplaceCanceller_ = defaultInplaceCanceller;
	if (!canceller) canceller_ = defaultCanceller;
	if (!observerAndReducer) observeAndReducer_ = defaultObserveReducer;
	if (!inplaceDamper) inplaceDamper_ = defaultInplaceWeakDamper;
} // Joint non-linear Gaussian constructor

CanonicalGaussianMixture::~CanonicalGaussianMixture() {} // Default Destructor

void CanonicalGaussianMixture::pruneComponents() {} // pruneComponents()

void CanonicalGaussianMixture::mergeComponents() {} // mergeComponents()

inline void CanonicalGaussianMixture::inplaceNormalize(FactorOperator* procPtr) {} // inplaceNormalize()

inline uniqptr<Factor> CanonicalGaussianMixture::normalize(FactorOperator* procPtr) const {
	return uniqptr<Factor> (new GaussCanonical());
} // normalize()

inline void CanonicalGaussianMixture::inplaceAbsorb(const Factor* rhsPtr, FactorOperator* procPtr) {
	if (procPtr) dynamicInplaceApply(procPtr, this, rhsPtr);
	else dynamicInplaceApply(inplaceAbsorber_.get(), this, rhsPtr);
} // inplaceAbsorb()

inline uniqptr<Factor> CanonicalGaussianMixture::absorb(const Factor* rhsPtr, FactorOperator* procPtr) const {
	return uniqptr<Factor> (new GaussCanonical());
} // absorb()

inline void CanonicalGaussianMixture::inplaceCancel(const Factor* rhsPtr, FactorOperator* procPtr) {} // inplaceCancel()

inline uniqptr<Factor> CanonicalGaussianMixture::cancel(const Factor* rhsPtr, FactorOperator* procPtr) const {
	return uniqptr<Factor> (new GaussCanonical());
} // cancel()

inline uniqptr<Factor> CanonicalGaussianMixture::marginalize(const emdw::RVIds& variablesToKeep, 
		bool presorted, FactorOperator* procPtr) const {
	return uniqptr<Factor> (new GaussCanonical());
} // marginalize()

inline uniqptr<Factor> CanonicalGaussianMixture::observeAndReduce( const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted, FactorOperator* procPtr) const {
	return uniqptr<Factor> (new GaussCanonical());
} // observeAndReduce()

double CanonicalGaussianMixture::inplaceDampen(const Factor* oldMsg, double df, FactorOperator* procPtr) {
	return 0.0;
} // inplaceDampen()

CanonicalGaussianMixture* CanonicalGaussianMixture::copy(const emdw::RVIds& newVars, bool presorted) const {
	return new CanonicalGaussianMixture(*this);
} // copy()

CanonicalGaussianMixture* CanonicalGaussianMixture::vacuousCopy(const emdw::RVIds& selectedVars, bool presorted) const {
	return new CanonicalGaussianMixture();
} // vacuousCopy()

bool CanonicalGaussianMixture::isEqual(const Factor* rhsPtr) const { return true; } // isEqual()

unsigned CanonicalGaussianMixture::noOfVars() const { return 0; } // noOfVars()

emdw::RVIds CanonicalGaussianMixture::getVars() const { return vars_; } // getVars()

emdw::RVIdType CanonicalGaussianMixture::getVar(unsigned varNo) const { return vars_[varNo]; } // getVar()

std::vector<rcptr<Factor>> CanonicalGaussianMixture::getComponents() const { return comps_; } //getComponents()

std::istream& CanonicalGaussianMixture::txtRead(std::istream& file) { return file; } // txtRead()

std::ostream& CanonicalGaussianMixture::txtWrite(std::ostream& file) const { return file; } // txtWrite()


//! Operators

const std::string& InplaceNormalizeCGM::isA() const {
	static const std::string CLASSNAME("InplaceNormalizeCGM");
	return CLASSNAME;
} // isA()

void InplaceNormalizeCGM::inplaceProcess(CanonicalGaussianMixture* lhsPtr) {
} // inplaceProcess()

const std::string& NormalizeCGM::isA() const {
	static const std::string CLASSNAME("NormalizeCGM");
	return CLASSNAME;
} // isA()

Factor* NormalizeCGM::process(const CanonicalGaussianMixture* lhsPtr) {
	return new CanonicalGaussianMixture();
} // process()

const std::string& InplaceAbsorbCGM::isA() const {
	static const std::string CLASSNAME("InplaceAbsorbCGM");
	return CLASSNAME;
} // isA()

void InplaceAbsorbCGM::inplaceProcess(CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {
	if (dynamic_cast<const CanonicalGaussianMixture*>(rhsFPtr) != NULL) {
		std::cout << "CanonicalGaussianMixture" << std::endl;
	} else {
		std::cout << "GaussCanonical" << std::endl;
	}
} // inplaceProcess()

const std::string& AbsorbCGM::isA() const {
	static const std::string CLASSNAME("AbsorbCGM");
	return CLASSNAME;
} // isA()

Factor* AbsorbCGM::process(const CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {
	return new CanonicalGaussianMixture();
} // process()

const std::string& InplaceCancelCGM::isA() const {
	static const std::string CLASSNAME("InplaceCancelCGM");
	return CLASSNAME;
} // isA()

void InplaceCancelCGM::inplaceProcess(CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {

} // inplaceCancel()

const std::string& CancelCGM::isA() const {
	static const std::string CLASSNAME("CancelCGM");
	return CLASSNAME;
} // isA()

Factor* CancelCGM::process(const CanonicalGaussianMixture* lhsPtr, const Factor* rhsFPtr) {
	return new CanonicalGaussianMixture();
} // process()

const std::string& MarginalizeCGM::isA() const {
	static const std::string CLASSNAME("MarginalizeCGM");
	return CLASSNAME;
} // isA()

Factor* MarginalizeCGM::process(const CanonicalGaussianMixture* lhsPtr, const emdw::RVIds& variablesToKeep,
		bool presorted) {
	return new CanonicalGaussianMixture();
} // process()

const std::string& ObserveAndReduceCGM::isA() const {
	static const std::string CLASSNAME("ObserveAndReduceCGM");
	return CLASSNAME;
} // isA()

Factor* ObserveAndReduceCGM::process(const CanonicalGaussianMixture* lhsPtr, const emdw::RVIds& variables,
		const emdw::RVVals& assignedVals, bool presorted) {
	return new CanonicalGaussianMixture();
} // process()

const std::string& InplaceWeakDampingCGM::isA() const {
	static const std::string CLASSNAME("InplaceWeakDampingCGM");
	return CLASSNAME;
} // isA()

double InplaceWeakDampingCGM::inplaceProcess(const CanonicalGaussianMixture* lhsPtr, const Factor* rhsPtr, double df) {
	return 0.0;
} // inplaceProcess()
