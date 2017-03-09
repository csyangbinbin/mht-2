/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for a Simplified Conditional Gaussian (Mixture) implementation.
 * See notes above class declaration.
 *************************************************************************/
#ifndef CONDITIONALGAUSSIAN_HPP  
#define CONDITIONALGAUSSIAN_HPP

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
class ConditionalGaussian;

/**
 * @brief Inplace normalization operator.
 */
class InplaceNormalizeCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		void inplaceProcess(ConditionalGaussian* lhsPtr);
}; // InplaceNormalizeCG

/**
 * @brief Normalization operator.
 */
class NormalizeCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const ConditionalGaussian* lhsPtr);
}; // NormalizeCG

/**
 * @brief Inplace absorbtion operator.
 */
class InplaceAbsorbCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		void inplaceProcess(ConditionalGaussian* lhsPtr,
				const Factor* rhsFPtr);
}; // InplaceAbsorbCG

/**
 * @brief Absorbtion operator.
 */
class AbsorbCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const ConditionalGaussian* lhsPtr,
				const Factor* rhsFPtr);
}; // AbsorbCG

/**
 * @brief Inplace cancellation operator.
 */
class InplaceCancelCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		void inplaceProcess(ConditionalGaussian* lhsPtr,
				const Factor* rhsFPtr);
}; // InplaceCancelCG

/**
 * @brief Cancellation operator.
 */
class CancelCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const ConditionalGaussian* lstPtr,
				const Factor* rhsFPtr);
}; // CancelCG

/**
 * @brief Marginalization operator.
 */
class MarginalizeCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const ConditionalGaussian* lhsPtr,
				const emdw::RVIds& variablesToKeep,
				bool presorted = false);
}; // MarginalizeCG

/**
 * @brief Observation and factor reduction operator.
 */
class ObserveAndReduceCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		Factor* process(const ConditionalGaussian* lhsPtr,
				const emdw::RVIds& variables,
				const emdw::RVVals& assignedVals,
				bool presorted = false);
}; // ObserveAndReduceCG

/**
 * @brief Inplace weak damping operator.
 */
class InplaceWeakDampingCG : public Operator1<ConditionalGaussian> {
	public:
		const std::string& isA() const;
		double inplaceProcess(const ConditionalGaussian* lhsPtr,
				const Factor* rhsPtr,
				double df);
}; // InplaceWeakDampingCG

/**
 * @brief Simplified Conditional Gaussian (Mixture)
 *
 * This is a simplified version of the CLG class.
 * The Gaussians are only allowed to be conditioned
 * on a single discrete RV. The prior over the
 * discrete varaible is assumed to absorbed into the
 * factor.
 *
 * This can be instantiated with any map of a 
 * discrete RV's domain to a Factor, but MarginalizeCG
 * requires either a GaussCanonical or CanonicalGaussian
 * Mixture.
 *
 * Not all required Factor methods are properly implemented, those
 * which aren't are noted in the documentation. inplaceWeakDamping,
 * txtRead and txtWrite are not implemented.
 *
 * To play it safe, there is a ridiculous amount of
 * copying of the conditional Factors it is probably
 * unnecessary.
 *
 * Never use in a ClusterGraph.
 *
 * @author SCJ Robertson
 * @since 10/02/17
 */
class ConditionalGaussian : public Factor {

	friend class InplaceNormalizeCG;
	friend class NormalizeCG;
	friend class InplaceAbsorbCG;
	friend class AbsorbCG;
	friend class InplaceCancelCG;
	friend class CancelCG;
	friend class MarginalizeCG;
	friend class ObserveAndReduceCG;
	friend class InplaceWeakDampingCG;

	public:
		/** 
		 * @brief Default vacuous constructor.
		 *  
		 * Default constrcutor.
		 *
		 * A list of Factor operators, if it equals zero it will be set to a default
		 * operator.
		 */
		ConditionalGaussian (
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
		ConditionalGaussian (
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

		ConditionalGaussian(const ConditionalGaussian& st) = default;
		
		ConditionalGaussian(ConditionalGaussian&& st) = default;

		/**
		 * @brief Default destructor.
		 */
		virtual ~ConditionalGaussian();

	public:
		ConditionalGaussian& operator=(const ConditionalGaussian& d) = default;

		ConditionalGaussian& operator=(ConditionalGaussian&& d) = default;


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
		virtual ConditionalGaussian* copy(const emdw::RVIds& newVars = {}, bool presorted = false ) const;

		/**
		 * @brief Vacuous copy
		 *
		 * Vacuous copy of the factor, possibly onto new or reduced scope.
		 */
		virtual ConditionalGaussian* vacuousCopy(const emdw::RVIds& selectedVars = {}, bool presorted = false) const;

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

		/**
		 * @brief Retuns the continuousVars
		 */
		emdw::RVIds getContinuousVars() const;

	public:
		/**
		 * @brief Read information from an input stream.
		 *	
		 * TODO: Implement this!
		 */
		virtual std::istream& txtRead(std::istream& file);

		/**
		 * @brief Write information to an output stream.
		 *
		 * TODO: Implement this!
		 */
		virtual std::ostream& txtWrite(std::ostream& file) const;

	// Data Members
	private:
		// Scope and components
		emdw::RVIds vars_;
		mutable std::map<unsigned, bool> isContinuous_; // Waste of space, but convenient

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

}; // ConditionalGaussian

#endif // CONDITIONALGAUSSIAN_HPP
