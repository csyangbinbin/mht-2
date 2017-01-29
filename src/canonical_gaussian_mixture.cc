/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for the Canonical Gaussian Mixture implementation.
 *************************************************************************/
#include <vector>
#include <iostream>
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
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& obseverAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {
	inplaceNormalizer_ = defaultInplaceNormalizer;
	normalizer_ = defaultNormalizer;
	inplaceAbsorber_ = defaultInplaceAbsorber;
	absorber_ = defaultAbsorber;
	inplaceCanceller_ = defaultInplaceCanceller;
	canceller_ = defaultCanceller;
	observeAndReducer_ = defaultObserveReducer;
	inplaceDamper_ = defaultInplaceWeakDamper;
} // Default Constructor

CanonicalGaussianMixture::CanonicalGaussianMixture(
		const emdw::RVIds& vars,
		const std::vector<double>& weights,
		const std::vector<ColVector<double>>& means,
		const std::vector<Matrix<double>>& covs,
		bool presorted,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& obseverAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {
	inplaceNormalizer_ = defaultInplaceNormalizer;
	normalizer_ = defaultNormalizer;
	inplaceAbsorber_ = defaultInplaceAbsorber;
	absorber_ = defaultAbsorber;
	inplaceCanceller_ = defaultInplaceCanceller;
	canceller_ = defaultCanceller;
	observeAndReducer_ = defaultObserveReducer;
	inplaceDamper_ = defaultInplaceWeakDamper;
} // Covariance constructor


CanonicalGaussianMixture::CanonicalGaussianMixture(
		const emdw::RVIds& vars,
		const std::vector<Matrix<double>>& prec,
		const std::vector<ColVector<double>>& info,
		const std::vector<double>& g,
		bool presorted,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& obseverAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {
	inplaceNormalizer_ = defaultInplaceNormalizer;
	normalizer_ = defaultNormalizer;
	inplaceAbsorber_ = defaultInplaceAbsorber;
	absorber_ = defaultAbsorber;
	inplaceCanceller_ = defaultInplaceCanceller;
	canceller_ = defaultCanceller;
	observeAndReducer_ = defaultObserveReducer;
	inplaceDamper_ = defaultInplaceWeakDamper;
} // Canonical constructor

/*
CanonicalGaussianMixture::CanonicalGaussianMixture(
		const emdw::RVIds& vars,
		const std::vector<rcptr<Factor>>& components,
		bool presorted = false,
		const rcptr<FactorOperator>& inplaceNormalizer,
		const rcptr<FactorOperator>& normalizer,
		const rcptr<FactorOperator>& inplaceAbsorber,
		const rcptr<FactorOperator>& absorber,
		const rcptr<FactorOperator>& inplaceCanceller,
		const rcptr<FactorOperator>& canceller,
		const rcptr<FactorOperator>& obseverAndReducer,
		const rcptr<FactorOperator>& inplaceDamper
		) {
	inplaceNormalizer_ = defaultInplaceNormalizer;
	normalizer_ = defaultNormalizer;
	inplaceAbsorber_ = defaultInplaceAbsorber;
	absorber_ = defaultAbsorber;
	inplaceCanceller_ = defaultInplaceCanceller;
	canceller_ = defaultCanceller;
	observeAndReducer_ = defaultObserveReducer;
	inplaceDamper_ = defaultInplaceWeakDamper;
} // Component constructor
*/

CanonicalGaussianMixture::~CanonicalGaussianMixture() {} // Default Destructor

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
	return new CanonicalGaussianMixture();
} // copy()

CanonicalGaussianMixture* CanonicalGaussianMixture::vacuousCopy(const emdw::RVIds& selectedVars, bool presorted) const {
	return new CanonicalGaussianMixture();
} // vacuousCopy()

bool CanonicalGaussianMixture::isEqual(const Factor* rhsPtr) const { return true; } // isEqual()

unsigned CanonicalGaussianMixture::noOfVars() const { return 0; } // noOfVars()

emdw::RVIds CanonicalGaussianMixture::getVars() const { return vars_; } // getVars()

emdw::RVIdType CanonicalGaussianMixture::getVar(unsigned varNo) const { return vars_[varNo]; } // getVar()

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
