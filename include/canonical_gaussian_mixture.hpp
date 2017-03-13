/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for the Canonical Gaussian Mixture implementation.
 * See the notes abve the class declaration.
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
 * @brief Prune the Gaussian mixture.
 *
 * Discard all components constributing insignificant mass. Also
 * discards anything with infinite mass, which is terrible
 * practice.
 *
 * @param components A vector of GaussCanonical factors representing
 * a GM.
 *
 * @param threshold The weight threshold.
 *
 * @return A reduced vector of GaussCanonical factors representing a GM.
 */
std::vector<rcptr<Factor>> pruneComponents( const std::vector<rcptr<Factor>>& components, const double threshold);

/**
 * @brief Merge all closely space components.
 *
 * Merge all components which are close to one another. This is
 * an incredibly expensive procedure, but there aren't any cheaper
 * methods of reducing a mixture.
 *
 * Typically used after pruneComponents. This is a
 * terrible implementation, but it seems to work. Some one
 * else check it?
 *
 * @param components A vector of GaussCanonical factors representing
 * a GM.
 *
 * @param maxComps The maximum of components allowed in a mixture.
 *
 * @param threshold The weight threshold.
 *
 * @param unionDistance The Mahalanobis distance between components
 *I
 * @return A reduced vector of GaussCanonical factors representing a GM.
 */
std::vector<rcptr<Factor>> mergeComponents( const std::vector<rcptr<Factor>>& components, const unsigned maxComps, 
		const double threshold, const double unionDistance);

/**
 * @brief Match a Gaussian Mixture's moments with a single Gaussian.
 *
 * Match a Gaussian mixture moments with a single GaussCanonical
 * factor.
 *
 * @param components A vector of GaussCanonical factors representing
 * a GM.
 *
 * @return A single GaussCanonical factor with moments matching the
 * GM
 */
uniqptr<Factor> mProject( const std::vector<rcptr<Factor>>& components);

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
 *
 * If the product has far too many components, pruneComponents
 * then mergeComponents is called.
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
 * for a vector of GaussCanonical factors, the weight of each
 * GM component is hidden in GaussCanonical's g component.
 *
 * For now, the underlying GaussCanonical factors make use of
 * their default operators. The mixture is limited in size
 * to some maximum number of components, after which insignificant 
 * components are pruned and the closely spaced components 
 * merged. If the mixture is still too large after pruning and merging, only
 * the N^{th} largest components are selected. See inplaceAbsorb,
 * pruneComponents and mergeComponents.
 *
 * Division of two mixtures will be approximated by first moment matching the
 * divisor and dividing each of the dividends components by the aprroxiamte Guassian.
 * Negative covariances can (and will) occur during division. Currently pruneComponents
 * throws out components with negative covariances (as they have infinite
 * mass and tend to break things.)
 *
 * Not all required Factor methods are properly implemented, those
 * which aren't are noted in the documentation. inplaceWeakDamping is
 * is not implemented.
 *
 * The vacuous constructor is alright in CanonicalGaussianMixture,
 * but GaussCanonical's vacuous constructor currently does 
 * initialize the variables properly and can lead to strange errors.
 * This probably because no one has ever used the vacuous constructor before.
 *
 * The implementation is fairly lax; it doesn't bother checking dimensional
 * consistency and the like. To play it safe, there is a ridiculous amount of
 * copying of the GM components which is expensive and probably unnecessary.
 *
 * Never use in a ClusterGraph.
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
	friend class MarginalizeCGM;
	friend class ObserveAndReduceCGM;
	friend class InplaceWeakDampingCGM;

	public:
		/** 
		 * @brief Default vacuous constructor.
		 * 
		 * Creates a single vacuous component Gaussian mixture.
		 *
		 * This currently is a bit dodgy from GaussCanonical's
		 * side, as its vacuous constructor doesn't initialise
		 * the g, h, K parameters as zero so weird errors may 
		 * result.
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
				const unsigned maxComponents = 3,
				const double threshold = 1e-6,
				const double unionDistance = 10,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& marginalizer = 0,
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
				const unsigned maxComponents = 3,
				const double threshold = 1e-6,
				const double unionDistance = 10,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& marginalizer = 0,
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
				const unsigned maxComponents = 3,
				const double threshold = 1e-6,
				const double unionDistance = 10,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& marginalizer = 0,
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
				const unsigned maxComponents = 3,
				const double threshold = 1e-6,
				const double unionDistance = 10,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& marginalizer = 0,
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
				const rcptr<Factor>& xFPtr,
				const Matrix<double>& A,
				const emdw::RVIds& newVars,
				const Matrix<double>& Q,
				bool presorted = false,
				const unsigned maxComponents = 3,
				const double threshold = 1e-6,
				const double unionDistance = 10,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& marginalizer = 0,
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
				const rcptr<Factor>& xFPtr,
				const rcptr<V2VTransform>& transform,
				const emdw::RVIds& newVars,
				const Matrix<double>& Q,
				bool presorted = false,
				const unsigned maxComponents = 3,
				const double threshold = 1e-6,
				const double unionDistance = 10,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& marginalizer = 0,
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
		virtual unsigned configure(unsigned key = 0);

		/**
		 * @brief Class specific configuration
		 *
		 * Reconfigures the exsiting clas, replacing old members with new 
		 * ones. It may seem pointless, but it greatly simplifies the 
		 * inplaceProcesses.
		 *list
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
		unsigned classSpecificConfigure(
				const emdw::RVIds& vars,
				const std::vector<rcptr<Factor>>& components,
				bool presorted = false,
				const unsigned maxComponents = 3,
				const double threshold = 1e-6,
				const double unionDistance = 10,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& marginalizer = 0,
				const rcptr<FactorOperator>& observerAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);
	public:
		/** 
		 * @brief Inplace normalization.
		 *
		 * Inplace normalization. Vacuous Gaussian distributions
		 * are unnormalizable, if this mixture contains a vacuous
		 * component normalization will fail.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 */
		inline void inplaceNormalize(FactorOperator* procPtr = 0);
		
		/**
		 * @brief Normalization.
		 *
		 * Normalization. Vacuous Gaussian distributions
		 * are unnormalizable, if this mixture contains a vacuous
		 * component normalization will fail.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 *
		 * @return A unique pointer to a normalized Gaussian mixture.
		 */
		inline uniqptr<Factor> normalize(FactorOperator* procPtr = 0) const;

		/**
		 * @brief Inplace multiplication.
		 *
		 * Inplace multiplication. The multiplier is allowed to be a CanonicalGaussianMixture
		 * or a GaussCanonical.
		 *
		 * @param rhsPtr The divisor. Can be CanonicalGaussianMixture or
		 * GaussCanonical.
		 */
		inline void inplaceAbsorb(const Factor* rhsPtr, FactorOperator* procPtr = 0);
		
		/**
		 * @brief Multiplication.
		 *
		 * Multiplication. The multiplier is allowed to be a CanonicalGaussianMixture
		 * or a GaussCanonical.
		 *
		 * @param rhsPtr The divisor. Can be CanonicalGaussianMixture or
		 * GaussCanonical.
		 *
		 * @return A uniqptr to the product Factor.
		 */
		inline uniqptr<Factor> absorb(const Factor* rhsPtr, FactorOperator* procPtr = 0) const;

		/**
		 * @brief Inplace division.
		 *
		 * Inplace division. The divisor is allowed to be a CanonicalGaussianMixture
		 * or a GaussCanonical. If it is a mixture, it will first be moment matched
		 * before division.
		 *
		 * @param rhsPtr The divisor. Can be CanonicalGaussianMixture or
		 * GaussCanonical.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 */
		inline void inplaceCancel(const Factor* rhsPtr, FactorOperator* procPtr = 0);
		
		/**
		 * @brief Division.
		 *
		 * Division. The divisor is allowed to be a CanonicalGaussianMixture
		 * or a GaussCanonical. If it is a mixture, it will first be moment matched
		 * before division.
		 *
		 * @param rhsPtr The divisor. Can be CanonicalGaussianMixture or
		 * GaussCanonical.
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 * process.
		 *
		 * @return A unique pointer to the quotient Factor.
		 */
		inline uniqptr<Factor> cancel(const Factor* rhsPtr, FactorOperator* procPtr = 0) const;
		
		/**
		 * @brief Marginalization.
		 *
		 * Marginalize out the given variables.
		 *
		 * @param variablesToKeep The variables which will not be marginalized out.
		 *
		 * @param presorted Is variablesToKeep sorted already?
		 *
		 * @param procPtr A unique pointer to some user defined FactorOperator
		 * process.
		 *
		 * @return A unique pointer to the scoped reduced Factor.
		 */
		inline uniqptr<Factor> marginalize(const emdw::RVIds& variablesToKeep, 
				bool presorted = false, FactorOperator* procPtr = 0) const;

		/**
		 * @brief Observe and Reduce
		 *
		 * Introduce evidence and collapse the factor.
		 *
		 * @param variables The variables which have been opbsevered in some 
		 * given state.
		 *
		 * @param assignedVals The values (state) of the given variables.
		 *
		 * @param presorted Are the given variables already sorted?
		 *
		 * @param procPtr A pointer to some user defined FactorOperator
		 *
		 * @return A unique pointer to resulting Factor.
		 */
		virtual uniqptr<Factor> observeAndReduce( const emdw::RVIds& variables,
				const emdw::RVVals& assignedVals, bool presorted = false,
				FactorOperator* procPtr = 0) const;

		/**
		 * @brief Inplace dampening.
		 *
		 * TODO: Complete this!!!
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
		 * @brief Moment match the mixture with a single Gaussian.
		 *
		 * Determine a single GaussCanonical which matches the mixture's
		 * first two moments.
		 *
		 * @return A unique pointer to a single GaussCanonical Factor
		 * with matching moments.
		 */
		uniqptr<Factor> momentMatch() const;

		/**
		 * @brief Prunes insignificant components and
		 * merges closely spaced components.
		 *
		 * Prunes insignificant components and merges
		 * closely spaced components. For practical
		 * purposes this must be called explicitly.
		 */
		void pruneAndMerge();

	public:
		/**
		 * @brief Adjust the mass of each component in the GM.
		 *
		 * Adjust the mass of each component in the GM. Essentially
		 * multiply the mixture by a constant.
		 *
		 * @param mass A constant shift of weight.
		 */
		void adjustMass(const double mass);

	public:
		/**
		 * @brief Return Gaussian mixture components.
		 */
		std::vector<rcptr<Factor>> getComponents() const;

		/**
		 * @brief Return number of Gaussian mixture components.
		 */
		double getNumberOfComponents() const;

		/**
		 * @brief Return the weights.
		 */
		std::vector<double> getWeights() const;

		/**
		 * @brief Return the means.
		 */
		std::vector<ColVector<double>> getMeans() const;

		/**
		 * @brief Return the covariance matrices.
		 */
		std::vector<Matrix<double>> getCovs() const;

		/**
		 * @brief Return the gs
		 */
		std::vector<double> getG() const;

		/**
		 * @brief Return the information vectors.
		 */
		std::vector<ColVector<double>> getH() const;

		/**
		 * @brief Return the precision matrices.
		 */
		std::vector<Matrix<double>> getK() const;

	// Data Members
	private:
		// Scope and components
		emdw::RVIds vars_;
		mutable std::vector<rcptr<Factor>> comps_;
		mutable unsigned N_;

		// Pruning and merging characteristics
		mutable unsigned maxComp_;
		mutable double threshold_;
		mutable double unionDistance_;

		// Operators
		rcptr<FactorOperator> inplaceNormalizer_;
		rcptr<FactorOperator> normalizer_;
		rcptr<FactorOperator> inplaceAbsorber_;
		rcptr<FactorOperator> absorber_;
		rcptr<FactorOperator> inplaceCanceller_;
		rcptr<FactorOperator> canceller_;
		rcptr<FactorOperator> marginalizer_;
		rcptr<FactorOperator> observeAndReducer_;
		rcptr<FactorOperator> inplaceDamper_;

}; // CanonicalGaussianMixture 

#endif // CANONCICALGAUSSIANMIXTURE_HPP
