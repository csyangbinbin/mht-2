/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for a Simplified Linear Gaussian implementation.
 *************************************************************************/
#ifndef LINEARGAUSSIAN_HPP
#define LINEARGAUSSIAN_HPP

#include <map>
#include "factor.hpp"
#include "factoroperator.hpp"
#include "emdw.hpp"
#include "anytype.hpp"
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "v2vtransform.hpp"

// Forward declaration.
class LinearGaussian;

/**
 * @brief Inplace normalization operator.
 */
class InplaceNormalizeLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		void inplaceProcess(LinearGaussian* lhsPtr);
}; // InplaceNormalizeLG

/**
 * @brief Normalization operator.
 */
class NormalizeLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const LinearGaussian* lhsPtr);
}; // NormalizeLG

/**
 * @brief Inplace absorbtion operator.
 */
class InplaceAbsorbLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		void inplaceProcess(LinearGaussian* lhsPtr,
				const Factor* rhsFPtr);
}; // InplaceAbsorbLG

/**
 * @brief Absorbtion operator.
 */
class AbsorbLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const LinearGaussian* lhsPtr,
				const Factor* rhsFPtr);
}; // AbsorbLG

/**
 * @brief Inplace cancellation operator.
 */
class InplaceCancelLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		void inplaceProcess(LinearGaussian* lhsPtr,
				const Factor* rhsFPtr);
}; // InplaceCancelLG

/**
 * @brief Cancellation operator.
 */
class CancelLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const LinearGaussian* lstPtr,
				const Factor* rhsFPtr);
}; // CancelLG

/**
 * @brief Marginalization operator.
 */
class MarginalizeLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const LinearGaussian* lhsPtr,
				const emdw::RVIds& variablesToKeep,
				bool presorted = false);
}; // MarginalizeLG

/**
 * @brief Observation and factor reduction operator.
 */
class ObserveAndReduceLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const LinearGaussian* lhsPtr,
				const emdw::RVIds& variables,
				const emdw::RVVals& assignedVals,
				bool presorted = false);
}; // ObserveAndReduceLG

/**
 * @brief Inplace weak damping operator.
 */
class InplaceWeakDampingLG : public Operator1<LinearGaussian> {
	public:
		const std::string& isA() const;
		double inplaceProcess(const LinearGaussian* lhsPtr,
				const Factor* rhsPtr,
				double df);
}; // InplaceWeakDampingLG

/**
 * @brief Simplified Conditional Linear Gaussian
 *
 * This is a simplified version of the CLG class.
 * The Gaussians are only allowed to be conditioned
 * on a single discrete RV. The prior over the
 * discrete varaible is assumed to abosrbed into the
 * factor.
 *
 * The CLG class in the emdw easily deal with 
 * GMM, so I recreated a simple class I understand.
 *
 * @author SCJ Robertson
 * @since 10/02/17
 */
class LinearGaussian : public Factor {

	friend class InplaceNormalizeLG;
	friend class NormalizeLG;
	friend class InplaceAbsorbLG;
	friend class AbsorbLG;
	friend class InplaceCancelLG;
	friend class CancelLG;
	friend class MarginalizeLG;
	friend class ObserveAndReduceLG;
	friend class InplaceWeakDampingLG;

	public:
		/** 
		 * @brief Default vacuous constructor.
		 *  
		 * Default constrcutor.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		LinearGaussian (
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
		 * @brief Class specific constructor
		 *
		 * Create a condtional linear Gaussian. The discrete RV is
		 * the prior over the conditioning variable. There must be
		 * a conditinuous fator assigned to every element in the discrete RV's
		 * domain. All continuous factors must have the same scope.
		 * 
		 * @param discreteRV A pointer to a distribution held over a
		 * single discrete random variable, it must be a DiscreteTable. 
		 * This is the prior held over the conditioning RV.
		 *
		 * @param conditionalList A map of the discrete variables domain
		 * to a continuous factor. Each factor in the list must have the same scope
		 * and the discrete RVs entire domain must map to a distribution.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		LinearGaussian (
				const rcptr<Factor>& discreteRV,
				const std::map<unsigned, rcptr<Factor>>& conditionalList,
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

		LinearGaussian(const LinearGaussian& st) = default;
		
		LinearGaussian(LinearGaussian&& st) = default;

		/**
		 * @brief Default destructor.
		 */
		virtual ~LinearGaussian();

	public:
		LinearGaussian& operator=(const LinearGaussian& d) = default;

		LinearGaussian& operator=(LinearGaussian&& d) = default;


	public:
		virtual unsigned configure(unsigned key = 0);

		/**
		 * @brief Class specific configuration
		 *
		 * Reconfigures the exsiting clas, replacing old members with new 
		 * ones. It may seem pointless, but it greatly simplifies the 
		 * inplaceProcesses.
		 *
		 * @param discreteRV A pointer to a distribution held over a
		 * single discrete random variable, it must be a DiscreteTable. 
		 * This is the prior held over the conditioning RV.
		 *
		 * @param conditionalList A map of the discrete variables domain
		 * to a continuous factor. Each factor in the list must have the same scope
		 * and the discrete RVs entire domain must map to a distribution.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		unsigned classSpecificConfigure(
				const rcptr<Factor>& discreteRV,
				const std::map<unsigned, rcptr<Factor>>& conditionalList,
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
		 * Division. The divisor is allowed to be a GaussCanonical or
		 * DiscreteTable.
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
		 * Division. The divisor is allowed to be a GaussCanonical or
		 * DiscreteTable.
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
		virtual LinearGaussian* copy(const emdw::RVIds& newVars = {}, bool presorted = false ) const;

		/**
		 * @brief Vacuous copy
		 *
		 * Vacuous copy of the factor, possibly onto new or reduced scope.
		 */
		virtual LinearGaussian* vacuousCopy(const emdw::RVIds& selectedVars = {}, bool presorted = false) const;

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

		/**
		 * @brief Return the prior over the discrete conditioning RV
		 */
		rcptr<Factor> getDiscretePrior() const;

		/**
		 * @brief Returns the continuous conditional distributions
		 */
		std::map<unsigned, rcptr<Factor>> getConditionalList() const;

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
		// Scope and components
		emdw::RVIds vars_;
		mutable std::map<unsigned, bool> isContinuous_; // Waste of space, but convinient

		// Factors
		rcptr<Factor> discreteRV_;
		mutable std::map<unsigned, rcptr<Factor>> conditionalList_;

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

}; // LinearGaussian

#endif // LINEARGAUSSIAN_HPP
