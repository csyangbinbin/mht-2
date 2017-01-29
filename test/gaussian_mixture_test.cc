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
		virtual void SetUp() {}

		virtual void TearDown() {}

};

TEST_F (CGMTest, DefaultConstructor) {
	rcptr<Factor> gm = uniqptr<CGM>(new CGM());
	EXPECT_EQ(0, 0);
}

TEST_F(CGMTest, InplaceAbsorber) {
	rcptr<Factor> gm_1 = uniqptr<CGM>(new CGM());
	rcptr<Factor> gm_2 = uniqptr<CGM>(new CGM());
	rcptr<Factor> gc = uniqptr<GaussCanonical>(new GaussCanonical());

	gm_1->inplaceAbsorb(gm_2);
	gm_1->inplaceAbsorb(gc);
}


