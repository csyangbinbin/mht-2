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
#include "linear_gaussian.hpp"
#include "utils.hpp"

class CLGTest : public testing::Test {

	protected:
		virtual void SetUp() {
		}

		virtual void TearDown() {
		}

	protected:
};
