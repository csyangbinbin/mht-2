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
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const emdw::RVIds& vars = {},
				bool presorted = false,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& obseverAndReducer = 0,
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
		 * must not be in logarithmic form.
		 *
		 * @param means A vector of the component's means.
		 *
		 * @param covs A vector of the components covariance matrices.
		 *
		 * @param presorted Set to true if vars is sorted according to their 
		 * integer values.
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
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& obseverAndReducer = 0,
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
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		CanonicalGaussianMixture(
				const emdw::RVIds& vars,
				const std::vector<Matrix<double>>& prec,
				const std::vector<ColVector<double>>& info,
				const std::vector<double>& g,
				bool presorted = false,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& obseverAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);


		/**
		 * @brief GaussCanonical vector constructor.
		 */
		/*
		CanonicalGaussianMixture(
				const emdw::RVIds& vars,
				const std::vector<rcptr<Factor>>& components,
				bool presorted = false,
				const rcptr<FactorOperator>& inplaceNormalizer = 0,
				const rcptr<FactorOperator>& normalizer = 0,
				const rcptr<FactorOperator>& inplaceAbsorber = 0,
				const rcptr<FactorOperator>& absorber = 0,
				const rcptr<FactorOperator>& inplaceCanceller = 0,
				const rcptr<FactorOperator>& canceller = 0,
				const rcptr<FactorOperator>& obseverAndReducer = 0,
				const rcptr<FactorOperator>& inplaceDamper = 0
				);
		*/


		CanonicalGaussianMixture(const CanonicalGaussianMixture& st) = default;
		CanonicalGaussianMixture(CanonicalGaussianMixture&& st) = default;

		/**
		 * @brief Default destructor.
		 */
		virtual ~CanonicalGaussianMixture();


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

	// Data Members
	private:
		emdw::RVIds vars_;
		std::vector<GaussCanonical> comps_;

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
