/*************************************************************************
 *  Compilation:  
 *  Execution:    
 *  Dependencies:
 *
 * Main app, runs small examples for now.
 *************************************************************************/

// patrec headers
#include "error.hpp"
#include "genmat.hpp"
#include "vecset.hpp"
#include "mean_cov.hpp"

// emdw headers
#include "emdw.hpp"
#include "anytype.hpp"
#include "clustergraph.hpp"
#include "gausscanonical.hpp"
#include "lbp_cg.hpp"
#include "matops.hpp"

// standard headers
#include <iostream>  // cout, endl, flush, cin, cerr
#include <cctype>  // toupper
#include <string>  // string
#include <vector>
#include <memory>
#include <set>
#include <map>
#include <limits>
#include <random>

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
	typedef GaussCanonical GCT1;
	std::vector<rcptr<Factor>> gaussians;
	
	enum{x_0, x_1, z_1, x_2, z_2, x_3, z_3};
	RVIds vars;

	vars = {x_0};
	gaussians.push_back( uniqptr<GCT1>( new GCT1(vars) ) );

	vars = {x_0, x_1};
	gaussians.push_back( uniqptr<GCT1>( new GCT1(vars) ) );

	vars = {x_1, z_1};
	gaussians.push_back( uniqptr<GCT1>( new GCT1(vars) ) );

	vars = {x_1, x_2};
	gaussians.push_back( uniqptr<GCT1>( new GCT1(vars) ) );

	vars = {x_2, z_2};
	gaussians.push_back( uniqptr<GCT1>( new GCT1(vars) ) );

	vars = {x_2, x_3};
	gaussians.push_back( uniqptr<GCT1>( new GCT1(vars) ) );

	vars = {x_3, z_3};
	gaussians.push_back( uniqptr<GCT1>( new GCT1(vars) ) );
	
	rcptr<ClusterGraph> cluster_graph;
	cluster_graph = uniqptr<ClusterGraph>(new ClusterGraph(gaussians));
	cluster_graph->exportToGraphViz("kalman");

	cout << "Not broken!" << endl;
}
