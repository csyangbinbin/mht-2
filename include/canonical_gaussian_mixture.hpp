/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for the Canonical Gaussian Mixture implementation.
 *************************************************************************/

#ifndef CANONICALGAUSSIANMIXTURE_HPP
#define CANONICALGAUSSIANMIXTURE_HPP

#include "factor.hpp"
#include "factoroperator.hpp"
#include "emdw.hpp"
#include "anytype.hpp"
#include "gausscanonical.hpp"
#include "v2vtransform.hpp"

// Forward declaration.
class CanonicalGaussianMixture;

/**
 * @brief Inplace normalization operator.
 */
class InplaceNormalizeCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		void inplaceProcess(CanonicalGaussianMixture* lhsPtr);
}; // InplaceNormalizeCGM

/**
 * @brief Normalization operator.
 */
class NormalizeCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		Factor* process(const CanonicalGaussianMixture* lhsPtr);
}; // NormalizeCGM

/**
 * @brief Inplace absorbtion operator.
 */
class InplaceAbsorbCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		void inplaceProcess(CanonicalGaussianMixture* lhsPtr,
				const Factor* rhsFPtr);
}; // InplaceAbsorbCGM

/**
 * @brief Absorbtion operator.
 */
class AbsorbCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		Factor* process(const CanonicalGaussianMixture* lhsPtr,
				const Factor* rhsFPtr);
}; // AbsorbCGM

/**
 * @brief Inplace cancellation operator.
 */
class InplaceCancelCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		void inplaceProcess(CanonicalGaussianMixture* lhsPtr,
				const Factor* rhsFPtr);
}; // InplaceCancelCGM

/**
 * @brief Cancellation operator.
 */
class CancelCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		Factor* process(const CanonicalGaussianMixture* lstPtr,
				const Factor* rhsFPtr);
}; // CancelCGM

/**
 * @brief Marginalization operator.
 */
class MarginalizeCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		Factor* process(const CanonicalGaussianMixture* lhsPtr,
				const emdw::RVIds& variablesToKeep,
				bool presorted = false);
}; // MarginalizeCGM

/**
 * @brief Observation and factor reduction operator.
 */
class ObserveAndReduceCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		Factor* process(const CanonicalGaussianMixture* lhsPtr,
				const emdw::RVIds& variables,
				const emdw::RVVals& assignedVals,
				bool presorted = false);
}; // ObserveAndReduceCGM

/**
 * @brief Inplace weak damping operator.
 */
class InplaceWeakDampingCGM : public Operator1<CanonicalGaussianMixture> {
	public:
		const std::string& isA() const;
		double inplaceProcess(const CanonicalGaussianMixture* lhsPtr,
				const Factor* rhsPtr,
				double df);
}; // InplaceWeakDampingCGM

/**
 * @brief Canonical Gaussian Mixture Model
 *
 * A rough Gaussian Mixture implementation. It is a wrapper type 
 * for a vector of GaussCanonical factors, it quietly pretends
 * NormedGaussCanonical doesn't exist. 
 *
 * Not all required Factor methods are properly implemented, those
 * which aren't are noted in the documentation.
 *
 * The implementation is fairly lax; it doesn't bother checking dimensional
 * consistency and the like.
 *
 * @author SCJ Robertson
 * @since 28/01/17
 */
class CanonicalGaussianMixture : public Factor {

	friend class InplaceNormalizeCGM;
	friend class NormalizeCGM;
	friend class InplaceAbsorbCGM;
	friend class AbsorbCGM;
	friend class InplaceCancelCGM;
	friend class CancelCGM;
	friend class ObserveAndReduceCGM;
	friend class InplaceWeakDampingCGM;
	friend class GaussCanonical;

	public:
		/** 
		 * @brief Default vacuous constructor.
		 * 
		 * Creates a single vacuous component Gaussian mixture.
		 * 
		 * @param vars Each variable in the PGM will be identified
		 * with a specific integer that indentifies it.
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
		 *
		 * @param maxComponents The maximum allowable number of components in the mixture.
		 *
		 * @param threshold The mimimum allowable mass a component is allowed to contribute.
		 *
		 * @param unionDistance The minimum Mahalanobis distance allowed between components.
		 * If the distance between their means is less than this threshold they merged into 
		 * one.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const emdw::RVIds& vars = {},
				bool presorted = false,
				const unsigned maxComponents = 100,
				const double threshold = 0.1,
				const double unionDistance = 5,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);

		/** 
		 * @brief Covariance form constructor.
		 * 
		 * Creates a Gaussian Mixture from the supplied vector of 
		 * means, covariances and weights.
		 *
		 * NOTE: The only difference between the covariance and
		 * canonical constructors is the order of the passed
		 * parameters.
		 *
		 * @param vars Each variable in the PGM will be identified
		 * with a specific integer that indentifies it.
		 *
		 * @param weights The weight of a mixture component. The weights
		 * must assume that Gaussian is otherwise normalised, they must given
		 * in linear form.
		 *
		 * @param means A vector of the component's means.
		 *
		 * @param covs A vector of the components covariance matrices.
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
		 *
		 * @param maxComponents The maximum allowable number of components in the mixture.
		 *
		 * @param threshold The mimimum allowable mass a component is allowed to contribute.
		 *
		 * @param unionDistance The minimum Mahalanobis distance allowed between components.
		 * If the distance between their means is less than this threshold they merged into 
		 * one.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const emdw::RVIds& vars,
				const std::vector<double>& weights,
				const std::vector<ColVector<double>>& means,
				const std::vector<Matrix<double>>& covs,
				bool presorted = false,
				const unsigned maxComponents = 100,
				const double threshold = 0.1,
				const double unionDistance = 5,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);

		/** 
		 * @brief Covariance form constructor.
		 * 
		 * Creates a Gaussian Mixture from the supplied vector of 
		 * information vectors, precision matrices and weights.
		 *
		 * NOTE: The only difference between the covariance and
		 * canonical constructors is the order of the passed
		 * parameters.
		 * 
		 * @param vars Each variable in the PGM will be identified
		 * with a specific integer that indentifies it.
		 *
		 * @param prec A vector of the components' precision matrices.
		 *
		 * @param info A vector of the components' information vectors.
		 *
		 * @param g A vector of the components' weights. The weights must be provided in
		 * logarithmic form.
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
		 *
		 * @param maxComponents The maximum allowable number of components in the mixture.
		 *
		 * @param threshold The mimimum allowable mass a component is allowed to contribute.
		 *
		 * @param unionDistance The minimum Mahalanobis distance allowed between components.
		 * If the distance between their means is less than this threshold they merged into 
		 * one.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const emdw::RVIds& vars,
				const std::vector<Matrix<double>>& prec,
				const std::vector<ColVector<double>>& info,
				const std::vector<double>& g,
				bool presorted = false,
				const unsigned maxComponents = 100,
				const double threshold = 0.1,
				const double unionDistance = 5,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);

		/**
		 * @brief GaussCanonical vector constructor.
		 *
		 * Creates a Gaussian mixture from the given vector of
		 * GaussCanonical factors.
		 *
		 * @param vars Each variable in the PGM will be identified
		 * with a specific integer that indentifies it.
		 *
		 * @param components A vector of GaussCanonical components, their
		 * variables must be the same as vars
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
		 *
		 * @param maxComponents The maximum allowable number of components in the mixture.
		 *
		 * @param threshold The mimimum allowable mass a component is allowed to contribute.
		 *
		 * @param unionDistance The minimum Mahalanobis distance allowed between components.
		 * If the distance between their means is less than this threshold they merged into 
		 * one.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const emdw::RVIds& vars,
				const std::vector<rcptr<Factor>>& components,
				bool presorted = false,
				const unsigned maxComponents = 100,
				const double threshold = 0.1,
				const double unionDistance = 5,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);

		/** 
		 * @brief Linear Gaussian constructor.
		 * 
		 * Creates a new Gaussian Mixture by shifting an exisiting
		 * Gaussian Mixture through a linear transfrom according to 
		 * Chapman-Kolmogrov equation.
		 * 
		 * @param xFPtr A pointer to an existing CanonicalGaussianMixture.
		 *
		 * @param A A matrix describing some appropriate linear transform.
		 *
		 * @param newVars The scope of the newly created GM.
		 *
		 * @param Q A noise covariance matrix.
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
		 *
		 * @param maxComponents The maximum allowable number of components in the mixture.
		 *
		 * @param threshold The mimimum allowable mass a component is allowed to contribute.
		 *
		 * @param unionDistance The minimum Mahalanobis distance allowed between components.
		 * If the distance between their means is less than this threshold they merged into 
		 * one.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const Factor* xFPtr,
				const Matrix<double>& A,
				const emdw::RVIds& newVars,
				const Matrix<double>& Q,
				bool presorted = false,
				const unsigned maxComponents = 100,
				const double threshold = 0.1,
				const double unionDistance = 5,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);

		/** 
		 * @brief Non-linear Gaussian constructor.
		 * 
		 * Creates a new Gaussian Mixture by shifting an exisiting
		 * Gaussian Mixture through a non-linear transfrom according to 
		 * Chapman-Kolmogrov equation using the Unscented Transform.
		 * 
		 * @param xFPtr A pointer to an existing CanonicalGaussianMixture.
		 *
		 * @param transfrom Any appropriate non-linear vector to vector
		 * transfrom.
		 *
		 * @param newVars The scope of the newly created GM.
		 *
		 * @param Q A noise covariance matrix.
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
		 *
		 * @param maxComponents The maximum allowable number of components in the mixture.
		 *
		 * @param threshold The mimimum allowable mass a component is allowed to contribute.
		 *
		 * @param unionDistance The minimum Mahalanobis distance allowed between components.
		 * If the distance between their means is less than this threshold they merged into 
		 * one.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const Factor* xFPtr,
				const V2VTransform& transform,
				const emdw::RVIds& newVars,
				const Matrix<double>& Q,
				bool presorted = false,
				const unsigned maxComponents = 100,
				const double threshold = 0.1,
				const double unionDistance = 5,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);

		/** 
		 * @brief Joint non-linear Gaussian constructor.
		 * 
		 * Creates a new joint Gaussian Mixture by shifting an exisiting
		 * Gaussian Mixture through a non-linear transfrom according to 
		 * Chapman-Kolmogrov equation using the Unscented Transform.
		 * 
		 * @param x1FPtr A pointer to an existing CanonicalGaussianMixture.
		 *
		 * @param x2FPtr A pointer to an existing CanonicalGaussianMixture.
		 *
		 * @param transfrom Any appropriate non-linear vector to vector
		 * transfrom.
		 *
		 * @param newVars The scope of the newly created GM.
		 *
		 * @param Q A noise covariance matrix.
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
		 *
		 * @param maxComponents The maximum allowable number of components in the mixture.
		 *
		 * @param threshold The mimimum allowable mass a component is allowed to contribute.
		 *
		 * @param unionDistance The minimum Mahalanobis distance allowed between components.
		 * If the distance between their means is less than this threshold they merged into 
		 * one.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const Factor* x1FPtr,
				const Factor* x2FPtr,
				const V2VTransform& transform,
				const emdw::RVIds& newVars,
				const Matrix<double>& Q,
				bool presorted = false,
				const unsigned maxComponents = 100,
				const double threshold = 0.1,
				const double unionDistance = 5,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);

		CanonicalGaussianMixture(const CanonicalGaussianMixture& st) = default;
		CanonicalGaussianMixture(CanonicalGaussianMixture&& st) = default;

		/**
		 * @brief Default destructor.
		 */
		virtual ~CanonicalGaussianMixture();

	public:
		CanonicalGaussianMixture& operator=(const CanonicalGaussianMixture& d) = default;
		CanonicalGaussianMixture& operator=(CanonicalGaussianMixture&& d) = default;


	public:
		/** 
		 * @brief Inplace normalization.
		 */
		inline void inplaceNormalize(FactorOperator* procPtr = 0);
		
		/**
		 * @brief Normalization.
		 */
		inline uniqptr<Factor> normalize(FactorOperator* procPtr = 0) const;

		/**
		 * @brief Inplace absorbtion.
		 */
		inline void inplaceAbsorb(const Factor* rhsPtr, FactorOperator* procPtr = 0);
		
		/**
		 * @brief Absorbtion.
		 */
		inline uniqptr<Factor> absorb(const Factor* rhsPtr, FactorOperator* procPtr = 0) const;

		/**
		 * @brief Inplace division.
		 */
		inline void inplaceCancel(const Factor* rhsPtr, FactorOperator* procPtr = 0);
		
		/**
		 * @brief Division.
		 */
		inline uniqptr<Factor> cancel(const Factor* rhsPtr, FactorOperator* procPtr = 0) const;
		
		/**
		 * @brief Marginalization.
		 */
		inline uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep, 
				bool presorted = false, FactorOperator* procPtr = 0) const;

		/**
		 * @brief Observe and Reduce
		 *
		 * Introduce evidence and collapse the factor.
		 */
		virtual uniqptr<Factor> observeAndReduce( const emdw::RVIds& variables,
				const emdw::RVVals& assignedVals, bool presorted = false,
				FactorOperator* procPtr = 0) const;

		/**
		 * @brief Inplace dampening.
		 */
		virtual double inplaceDampen(const Factor* oldMsg, double df, FactorOperator* procPtr = 0);

	public:
		/** 
		 * @brief Copy factor.
		 *
		 * Copy factor, possibly onto new scope.
		 */
		virtual CanonicalGaussianMixture* copy(const emdw::RVIds& newVars = {}, bool presorted = false ) const;

		/**
		 * @brief Vacuous copy
		 *
		 * Vacuous copy of the factor, possibly onto new or reduced scope.
		 */
		virtual CanonicalGaussianMixture* vacuousCopy(const emdw::RVIds& selectedVars = {}, bool presorted = false) const;

		/**
		 * @brief Equality.
		 */
		bool isEqual(const Factor* rhsPtr) const;

		/**
		 * @brief Distance from a vacuous distribution.
		 *
		 * Factor distance from a vacuous distribution. Lets the base class
		 * take care of implementation.
		 */
		double distanceFromVacuous() const { return Factor::distanceFromVacuous(); }

		/**
		 * @brief Returns number of variables.
		 */
		virtual unsigned noOfVars() const;

		/**
		 * @brief Returns the variables' identity.
		 */
	virtual emdw::RVIds getVars() const;

		/**
		 * @brief Returns a specfic variable's identity.
		 */
		virtual emdw::RVIdType getVar(unsigned varNo) const;

	public:
		/**
		 * @brief Read information from an input stream.
		 */
		virtual std::istream& txtRead(std::istream& file);

		/**
		 * @brief Write information to an output stream.
		 */
		virtual std::ostream& txtWrite(std::ostream& file) const;

	public:
		/**
		 * @brief Prune the Gaussian mixture.
		 *
		 * Discard all components constributing insignificant mass.
		 */
		void pruneComponents();

		/**
		 * @brief Merge all closely space components.
		 *
		 * Merge all components which are close to one another.
		 */
		void mergeComponents();


	// Data Members
	private:
		// Scope and components
		emdw::RVIds vars_;
		std::vector<rcptr<Factor>> comps_;
		mutable unsigned N_;

		// Canonical representation
		// Not sure if these are actually worth keeping around, all info is already in comps_
		mutable std::vector<double> g_;
		mutable std::vector<ColVector<double>> h_;
		mutable std::vector<Matrix<double>> K_;

		// Pruning and merging characteristics
		mutable unsigned maxComp_;
		mutable double threshold_;
		mutable double unionDistance_;

		// Operators
		rcptr<FactorOperator> marginalizer_;
		rcptr<FactorOperator> inplaceNormalizer_;
		rcptr<FactorOperator> normalizer_;
		rcptr<FactorOperator> inplaceAbsorber_;
		rcptr<FactorOperator> absorber_;
		rcptr<FactorOperator> inplaceCanceller_;
		rcptr<FactorOperator> canceller_;
		rcptr<FactorOperator> observeAndReducer_;
		rcptr<FactorOperator> inplaceDamper_;

}; // CanonicalGaussianMixture 

#endif // CANONCICALGAUSSIANMIXTURE_HPP
