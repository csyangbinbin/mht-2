/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for canonical_gaussian_mixture.hpp.
 *************************************************************************/
#include <iostream>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "anytype.hpp"
#include "emdw.hpp"
#include "oddsandsods.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "transforms.hpp"
#include "utils.hpp"

using namespace emdw;
using namespace std;

class CGMTest : public testing::Test {

	protected:
		typedef CanonicalGaussianMixture CGM;

	protected:
		virtual void SetUp() {
			// Identity
			vars_ = emdw::RVIds(kDim_);
			for (unsigned i = 0; i < kDim_; i++) vars_[i] = i;
			
			// Weights
			w_ = std::vector<double>(kCompN_);
			for (unsigned i = 0; i < kCompN_; i++) w_[i] = 1.0/kCompN_;

			// Mean and Covariances
			mu_ = std::vector<ColVector<double>>(kCompN_);
			S_ = std::vector<Matrix<double>>(kCompN_);
			for (unsigned i = 0; i < kCompN_; i++) {
				mu_[i] = ColVector<double>(kDim_); mu_[i] *= 0;

				S_[i] = gLinear::zeros<double> (kDim_, kDim_);
				for (unsigned j = 0; j < kDim_; j++) (S_[i])(j, j) = 1;
			}

			// Canonical Form representation
			g_ = std::vector<double>(kCompN_);
			h_ = std::vector<ColVector<double>>(kCompN_);
			K_ = std::vector<Matrix<double>>(kCompN_);

			double detcov;
			int fail;
			for (unsigned i = 0; i < kCompN_; i++) {
				K_[i] = inv(S_[i], detcov, fail);
				if (fail) printf("Could not invert S_[%d] at line number %d in file %s\n",
						i, __LINE__, __FILE__);

				h_[i] = S_[i]*mu_[i];
				g_[i] = -0.5*( (K_[i]*mu_[i]).transpose() )*mu_[i] 
					- log( ( pow(2*M_PI, (1.0*vars_.size())/2) * pow(detcov, 0.5) ) / w_[i]);
			}

		}

		virtual void TearDown() {
			vars_.clear();

			w_.clear();
			mu_.clear();
			S_.clear();

			g_.clear();
			h_.clear();
			K_.clear();
		}

	protected:
		const unsigned kCompN_ = 3;
		const unsigned kDim_ = 3;
		emdw::RVIds vars_;

		std::vector<double> w_;
		std::vector<ColVector<double>> mu_;
		std::vector<Matrix<double>> S_;

		std::vector<double> g_;
		std::vector<ColVector<double>> h_;
		std::vector<Matrix<double>> K_;

		const unsigned kMaxComp_ = 100;
		const double kThreshold_ = 0.1;
		const double kUnionDistance_ = 5; 
};

TEST_F (CGMTest, DefaultConstructor) {
	rcptr<Factor> gm = uniqptr<CGM>(new CGM());
}

TEST_F (CGMTest, CovarianceConstructor) {
	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, w_, mu_, S_, false));
}

TEST_F (CGMTest, CanonicalConstructor) {
	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_, false));
}

TEST_F (CGMTest, ComponentConstructor) {
	std::vector<rcptr<Factor>> gm_vec = std::vector<rcptr<Factor>> (kCompN_);
	for (unsigned i = 0; i < kCompN_; i++) gm_vec[i] = uniqptr<GaussCanonical>(new GaussCanonical(vars_, K_[i], h_[i], g_[i]));

	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, gm_vec));
}

TEST_F(CGMTest, InplaceAbsorber) {
	rcptr<Factor> gm_1 = uniqptr<CGM>(new CGM());
	rcptr<Factor> gm_2 = uniqptr<CGM>(new CGM());
	rcptr<Factor> gc = uniqptr<GaussCanonical>(new GaussCanonical());

	gm_1->inplaceAbsorb(gm_2);
	gm_1->inplaceAbsorb(gc);
}


