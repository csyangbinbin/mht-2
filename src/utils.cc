/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * A library of useful helper functions which don't really fit in anywhere 
 * else. 
 *************************************************************************/
#include "utils.hpp"

unsigned addVariables (emdw::RVIds& globalVariables,
		emdw::RVIds& localVariables, 
		std::map<unsigned, emdw::RVIds>& map, 
		const unsigned L) {
	// Temporary variables
	unsigned M, N;
	emdw::RVIds temp; temp.clear(); temp.resize(L);

	// Add to the local variables
	M = localVariables.size(); localVariables.resize(M + 1);
	localVariables[M] = M;

	// Add to the global variables
	N = globalVariables.size(); globalVariables.resize(N + L);
	for (unsigned i = N; i < N+L; i++) {
		globalVariables[N] = i;
		temp[i - N] = i;
	}

	// Add to the map
	map[M].clear();
	map[M] = temp;

	return localVariables[M];
} // addVariables()
