/*************************************************************************
 *  Compilation:  
 *  Execution:    
 *  Dependencies:
 *
 * Main app, runs small examples for now.
 *************************************************************************/

// patrec headers
#include "sparsevec.hpp"

// emdw headers
#include "emdw.hpp"
#include "anytype.hpp"
#include "clustergraph.hpp"
#include "gausscanonical.hpp"
#include "oddsandsods.hpp"
#include "clg.hpp"
#include "lbp_cg.hpp"
#include "matops.hpp"

// standard headers
#include <iostream>  
#include <string>  
#include <map>
#include <vector>

#include "v2vtransform.hpp"

using namespace emdw;
using namespace std;

/**
 * Main app, runs small examples for now.
 *
 * @author SCJ Robertson
 * @since 0.9.1
 */
int main(int, char *argv[]) {

	typedef GaussCanonical GC;
	typedef ConditionalGauss CG;
	
	enum{a_0, a_1, a_2, x0_0, x0_1, z0_1, z1_1};

	std::vector<rcptr<Factor>> factors;
	std::map<std::vector<unsigned>, rcptr<Factor>> gaussian_map;

	//Association structure
	RVIds discrete_vars{a_1, a_2};
	RVIds cont_vars{x0_1, z0_1, z1_1};

	SparseVector<std::vector<unsigned>, double> assoc_prob;
	assoc_prob[RVIds{a_1}] = 0.5;
	assoc_prob[RVIds{a_2}] = 0.5;

	rcptr<Factor> likelihood_1 = uniqptr<GC>(new GC(RVIds{x0_1, z0_1}));
	rcptr<Factor> likelihood_2 = uniqptr<GC>(new GC(RVIds{x0_1, z1_1}));

	gaussian_map.clear();
	gaussian_map[RVIds{a_1}] = likelihood_1;
	gaussian_map[RVIds{a_2}] = likelihood_2;

	rcptr<Factor> markovian_1 = uniqptr<GC>(new GC(RVIds{x0_0, x0_1}));
	rcptr<Factor> association_1 = uniqptr<CG>(new CG(discrete_vars, assoc_prob, cont_vars, gaussian_map));

	factors.push_back(markovian_1);
	//factors.push_back(association_1);
	
	rcptr<ClusterGraph> cluster_graph = uniqptr<ClusterGraph> (new ClusterGraph(factors));

	cout << discrete_vars << endl;
}
