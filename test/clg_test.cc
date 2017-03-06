/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for linear_gaussian.hpp.
 *************************************************************************/
#include <iostream>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "anytype.hpp"
#include "emdw.hpp"
#include "oddsandsods.hpp"
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "linear_gaussian.hpp"
#include "utils.hpp"

class CLGTest : public testing::Test {

	protected:
		virtual void SetUp() {
			inplaceNormalizer_ = uniqptr<FactorOperator> (new DiscreteTable_InplaceNormalize<T>);
			normalizer_ = uniqptr<FactorOperator> (new DiscreteTable_Normalize<T>);
			marginalizer_ = uniqptr<FactorOperator> (new DiscreteTable_Marginalize<T>);
			inplaceCanceller_ = uniqptr<FactorOperator> (new DiscreteTable_InplaceCancel<T>);
			canceller_ = uniqptr<FactorOperator> (new DiscreteTable_Cancel<T>);

			a0Dom_ = uniqptr<DASS>(new DASS{0, 1});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{1}] = sparseProbs_[DASS{2}] = 1;

			discreteRV_ = uniqptr<Factor> (new DT(emdw::RVIds{a0}, {a0Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, marginalizer_, inplaceNormalizer_, normalizer_) );

			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{1}] = sparseProbs_[DASS{2}] = 2;
			discreteMsg_ = uniqptr<Factor> (new DT(emdw::RVIds{a0}, {a0Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, marginalizer_, inplaceNormalizer_, normalizer_) );

			// Create the conditional list
			mu_ = std::vector<ColVector<double>>(kCompN_);
			S_ = std::vector<Matrix<double>>(kCompN_);
			continuousFactors_ = std::vector<rcptr<Factor>>(kCompN_);

			for (unsigned i = 0; i < kCompN_; i++) {
				mu_[i] = ColVector<double>(kDim_); mu_[i] *= 0;

				S_[i] = gLinear::zeros<double> (kDim_, kDim_);
				for (unsigned j = 0; j < kDim_; j++) (S_[i])(j, j) = 1;

				continuousFactors_[i] = uniqptr<Factor>(new GaussCanonical(emdw::RVIds{x0, x1}, mu_[i], S_[i]));
				conditionalList_[T(i)] = uniqptr<Factor>(continuousFactors_[i]->copy());
			}

			// Create a GMM
			h_ = std::vector<ColVector<double>>(kCompN_);
			K_ = std::vector<Matrix<double>>(kCompN_);
			mixtureComponents_ = std::vector<rcptr<Factor>>(kCompN_);

			for (unsigned i = 0; i < kCompN_; i++) {
				h_[i] = ColVector<double>(kDim_); h_[i] *= 0;
				K_[i] = gLinear::zeros<double> (kDim_, kDim_);
				for (unsigned j = 0; j < kDim_; j++) {
					(h_[i])[j] = 2;
					(K_[i])(j, j) = 2;
				} 
				mixtureComponents_[i] = uniqptr<Factor>(new GaussCanonical(emdw::RVIds{x0, x1}, K_[i], h_[i]));
			}
			gm_ = uniqptr<Factor>(new CanonicalGaussianMixture(emdw::RVIds{x0, x1}, mixtureComponents_));
		}

		virtual void TearDown() {
			sparseProbs_.clear();
			mu_.clear();
			S_.clear();
			conditionalList_.clear();
		}


	protected:
		typedef unsigned short T;
		typedef DiscreteTable<T> DT;
		typedef std::vector<T> DASS;

	protected:
		// Vars
		enum{a0, x0, x1};

		// Discrete
		const double kFloor_ = 0.0;
		const double kMargin_ = 0.0;
		const double kDefProb_ = 0.0;

		rcptr<FactorOperator> inplaceNormalizer_;
		rcptr<FactorOperator> normalizer_;
		rcptr<FactorOperator> marginalizer_;
		rcptr<FactorOperator> inplaceCanceller_;
		rcptr<FactorOperator> canceller_;

		rcptr<DASS> a0Dom_;
		std::map<DASS, FProb> sparseProbs_;
		rcptr<Factor> discreteRV_;
		rcptr<Factor> discreteMsg_;

		// Continuous
		const unsigned kCompN_ = 3;
		const unsigned kDim_ = 2;

		std::vector<ColVector<double>> mu_;
		std::vector<Matrix<double>> S_;

		std::vector<ColVector<double>> h_;
		std::vector<Matrix<double>> K_;

		std::vector<rcptr<Factor>> continuousFactors_;
		std::vector<rcptr<Factor>> mixtureComponents_;
		rcptr<Factor> gm_;

		// Map
		std::map<unsigned, rcptr<Factor>> conditionalList_;
};

TEST_F (CLGTest, InitTest) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
}

TEST_F (CLGTest, InplaceNormalize) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	lg->inplaceNormalize();
}

TEST_F (CLGTest, InplaceAbsorbDiscrete) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	lg->inplaceAbsorb(discreteMsg_);

	rcptr<LinearGaussian> cast = std::dynamic_pointer_cast<LinearGaussian>(lg);
	rcptr<Factor> discretePrior = cast->getDiscretePrior();
}

TEST_F (CLGTest, InplaceAbsorbGaussian) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	rcptr<Factor> gc = uniqptr<Factor>( (mixtureComponents_[0])->copy() );
	lg->inplaceAbsorb( gc );

	rcptr<LinearGaussian> cast = std::dynamic_pointer_cast<LinearGaussian>(lg);
	std::map<unsigned, rcptr<Factor>> map = cast->getConditionalList();
}

TEST_F (CLGTest, InplaceAbsorbGM) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	lg->inplaceAbsorb(gm_);

	rcptr<LinearGaussian> cast = std::dynamic_pointer_cast<LinearGaussian>(lg);
	std::map<unsigned, rcptr<Factor>> map = cast->getConditionalList();
	rcptr<CanonicalGaussianMixture> cgm;

	for (auto& i : map) cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(i.second);
}

TEST_F (CLGTest, InplaceCancelDiscrete) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	lg->inplaceCancel(discreteMsg_);

	rcptr<LinearGaussian> cast = std::dynamic_pointer_cast<LinearGaussian>(lg);
	rcptr<Factor> discretePrior = cast->getDiscretePrior();
}

TEST_F (CLGTest, InplaceCancelGaussian) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	rcptr<Factor> gc = uniqptr<Factor>( (mixtureComponents_[0])->copy() );
	lg->inplaceCancel( gc );

	rcptr<LinearGaussian> cast = std::dynamic_pointer_cast<LinearGaussian>(lg);
	std::map<unsigned, rcptr<Factor>> map = cast->getConditionalList();
}

TEST_F (CLGTest, InplaceCancelGM) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	lg->inplaceCancel(gm_);

	rcptr<LinearGaussian> cast = std::dynamic_pointer_cast<LinearGaussian>(lg);
	std::map<unsigned, rcptr<Factor>> map = cast->getConditionalList();
}

TEST_F (CLGTest, MarginalizeDiscrete) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
	rcptr<Factor> gm = lg->marginalize(emdw::RVIds{x0, x1}); 

	rcptr<CanonicalGaussianMixture> cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(gm);
	std::vector<rcptr<Factor>> comps = cgm->getComponents();
}

TEST_F (CLGTest, ObserveAndReduceDiscrete) {
	rcptr<Factor> lg = uniqptr<Factor>(new LinearGaussian(discreteRV_, conditionalList_));
}
