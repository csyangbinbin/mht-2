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
			for (unsigned i = 0; i < kCompN_; i++) w_[i] = 1.0;
			w_[0] = 5; w_[kCompN_ - 1] = 0.05;

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

			// Linear motion model
			A_ = gLinear::zeros<double> (kDim_, kDim_);
			for (unsigned i = 0; i < kDim_; i++) A_(i, i) = 1;

			// Noise
			Q_ = gLinear::zeros<double> (kDim_, kDim_);
			for (unsigned i = 0; i < kDim_; i++) Q_(i, i) = 1;

			// MotionModel
			motion_model_ = uniqptr<MotionModel>(new MotionModel(kTimeStep_));
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
		const unsigned kDim_ = 6;
		emdw::RVIds vars_;

		std::vector<double> w_;
		std::vector<ColVector<double>> mu_;
		std::vector<Matrix<double>> S_;

		std::vector<double> g_;
		std::vector<ColVector<double>> h_;
		std::vector<Matrix<double>> K_;

		Matrix<double> A_;
		Matrix<double> Q_;

		rcptr<V2VTransform> motion_model_;
		rcptr<V2VTransform> sensor_model_;

		const unsigned kMaxComp_ = 100;
		const double kThreshold_ = 0.001;
		const double kUnionDistance_ = 5; 
		const double kTimeStep_ = 0.04;
};

TEST_F (CGMTest, DefaultConstructor) {
	rcptr<CGM> gm = uniqptr<CGM>(new CGM());
}

TEST_F (CGMTest, CovarianceConstructor) {
	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, w_, mu_, S_));
}

TEST_F (CGMTest, CanonicalConstructor) {
	rcptr<CGM> gm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
}

TEST_F (CGMTest, ComponentConstructor) {
	std::vector<rcptr<Factor>> gm_vec = std::vector<rcptr<Factor>> (kCompN_);
	for (unsigned i = 0; i < kCompN_; i++) gm_vec[i] = uniqptr<GaussCanonical>(new GaussCanonical(vars_, K_[i], h_[i], g_[i]));

	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, gm_vec));
}

TEST_F (CGMTest, LinearConstructor) {
	emdw::RVIds newVars = emdw::RVIds(kDim_);
	rcptr<Factor> oldGM = uniqptr<CGM>(new CGM(vars_, K_, h_, g_, false));

	for (unsigned i = 0; i < kDim_; i++) newVars[i] = kDim_ + i;
	rcptr<Factor> newGM = uniqptr<CGM>(new CGM(oldGM, A_, newVars, Q_));
}

TEST_F (CGMTest, NonlinearConstructor) {
	emdw::RVIds newVars = emdw::RVIds(kDim_);
	rcptr<Factor> oldGM = uniqptr<CGM>(new CGM(vars_, K_, h_, g_, false));

	for (unsigned i = 0; i < kDim_; i++) newVars[i] = kDim_ + i;
	rcptr<Factor> newGM = uniqptr<CGM>(new CGM(oldGM, motion_model_, newVars, Q_));
}

TEST_F (CGMTest, InplaceAbsorbGC) {
	emdw::RVIds newVars = emdw::RVIds(kDim_);
	for (unsigned i = 0; i < kDim_; i++) newVars[i] = kDim_ + i;

	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
	rcptr<Factor> gc = uniqptr<GaussCanonical>(new GaussCanonical(newVars, K_[0], h_[0], g_[0]));
	gm->inplaceAbsorb(gc);
}

TEST_F (CGMTest, InplaceAbsorbCGM) {
	emdw::RVIds newVars = emdw::RVIds(kDim_);
	for (unsigned i = 0; i < kDim_; i++) newVars[i] = kDim_ + i;

	rcptr<Factor> multiplicand = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
	rcptr<Factor> multiplier = uniqptr<CGM>(new CGM(newVars, K_, h_, g_));

	multiplicand->inplaceAbsorb(multiplier);
} 

TEST_F (CGMTest, AbsorbGC) {
	emdw::RVIds newVars = emdw::RVIds(kDim_);
	for (unsigned i = 0; i < kDim_; i++) newVars[i] = kDim_ + i;

	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
	rcptr<Factor> gc = uniqptr<GaussCanonical>(new GaussCanonical(newVars, K_[0], h_[0], g_[0]));
	rcptr<Factor> product = gm->absorb(gc);
}


TEST_F (CGMTest, AbsorbGCM) {
	emdw::RVIds newVars = emdw::RVIds(kDim_);
	for (unsigned i = 0; i < kDim_; i++) newVars[i] = kDim_ + i;

	rcptr<Factor> multiplicand = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
	rcptr<Factor> multiplier = uniqptr<CGM>(new CGM(newVars, K_, h_, g_));
	rcptr<Factor> product = multiplicand->absorb(multiplier);
}

TEST_F (CGMTest, InplaceNormalize) {
	rcptr<CGM> cgm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_, false, kMaxComp_, kThreshold_, kUnionDistance_));
	std::vector<rcptr<Factor>> cgmComps = cgm->getComponents();
	
	cgm->inplaceNormalize();
}

TEST_F (CGMTest, Marginalize) {
	rcptr<Factor> gm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
	rcptr<Factor> reduced = gm->marginalize(emdw::RVIds{1, 2});
}

TEST_F (CGMTest, PruneComponents) {
	rcptr<CGM> gm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
	gm->pruneComponents();
}

TEST_F (CGMTest, MergeComponents) {
	rcptr<CGM> gm = uniqptr<CGM>(new CGM(vars_, K_, h_, g_));
	gm->mergeComponents();
	
	std::vector<rcptr<Factor>> comps = gm->getComponents();
	//std::cout << *comps[0] << std::endl;
}
