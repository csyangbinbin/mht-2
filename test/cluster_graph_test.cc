/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Google Test fixture for graph.hpp and node.hpp.
 *************************************************************************/
#include <iostream>
#include <map>
#include "gtest/gtest.h"
#include "genvec.hpp"
#include "genmat.hpp"
#include "anytype.hpp"
#include "emdw.hpp"
#include "oddsandsods.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "node.hpp"
#include "transforms.hpp"
#include "utils.hpp"

using namespace emdw;
using namespace std;

class GraphTest : public testing::Test {
	protected:
		virtual void SetUp() {
		}

		virtual void TearDown() {
		}
};

TEST_F (GraphTest, InitTest) {
	std::map <rcptr<Factor>, emdw::RVIds> m;
	m.clear();

	rcptr<Factor> gc_key = uniqptr<GaussCanonical>(new GaussCanonical());
	emdw::RVIds val = {0, 1, 3};

	m[gc_key] = val;

	std::cout << m[gc_key] << std::endl;
}
