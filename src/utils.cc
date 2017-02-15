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
	unsigned M, N;

	M = localVariables.size();
	localVariables.push_back(M);

	N = globalVariables.size();
	for (unsigned i = 0; i < L; i++) globalVariables.push_back(N + i);

	map[M] = emdw::RVIds(&globalVariables[N+L], &globalVariables[N+L] + L);

	return localVariables[M];
} // addVariables()
