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

			a0Dom_ = uniqptr<DASS>(new DASS{0, 1});
			sparseProbs_.clear();
			sparseProbs_[DASS{0}] = sparseProbs_[DASS{1}] = sparseProbs_[DASS{2}] = 1;

			discreteRV_ = uniqptr<DT> (new DT(emdw::RVIds{a0}, {a0Dom_}, kDefProb_, 
					sparseProbs_, kMargin_, kFloor_, false, marginalizer_, inplaceNormalizer_, normalizer_) );

			mu_ = std::vector<ColVector<double>>(kCompN_);
			S_ = std::vector<Matrix<double>>(kCompN_);
			continuousFactors_ = std::vector<rcptr<Factor>>(kCompN_);

			for (unsigned i = 0; i < kCompN_; i++) {
				mu_[i] = ColVector<double>(kDim_); mu_[i] *= 0;

				S_[i] = gLinear::zeros<double> (kDim_, kDim_);
				for (unsigned j = 0; j < kDim_; j++) (S_[i])(j, j) = 1;

				continuousFactors_[i] = uniqptr<Factor>(new GaussCanonical(emdw::RVIds{x0, x1}, mu_[i], S_[i]));
				conditionalList_[T(i)] = continuousFactors_[i];
			}
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

		rcptr<DASS> a0Dom_;
		std::map<DASS, FProb> sparseProbs_;
		rcptr<Factor> discreteRV_;

		// Continuous
		const unsigned kCompN_ = 3;
		const unsigned kDim_ = 2;

		std::vector<ColVector<double>> mu_;
		std::vector<Matrix<double>> S_;

		std::vector<rcptr<Factor>> continuousFactors_;

		// Map
		std::map<unsigned, rcptr<Factor>> conditionalList_;
};

TEST_F (CLGTest, InitTest) {
	rcptr<LinearGaussian> lg = uniqptr<LinearGaussian>(new LinearGaussian(discreteRV_, conditionalList_));
}

TEST_F (CLGTest, InplaceNormalize) {
	rcptr<LinearGaussian> lg = uniqptr<LinearGaussian>(new LinearGaussian(discreteRV_, conditionalList_));

	//lg->classSpecificConfigure(discreteRV_, conditionalList_);
	
	lg->inplaceNormalize();

	//std::cout << lg->getDiscretePrior() << std::endl;
}
