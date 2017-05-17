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

double calculateEvidence(const unsigned K, 
		const unsigned N,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	
	double odds = 0;

	for (unsigned i = K; i <= N; i++) {
		unsigned M = stateNodes[i].size();
		
		// Determine the log-odds - including vacuous sponge.
		for (unsigned j = 0; j < M; j++) {
			if (stateNodes[N][j] == 0) continue;
			double mass = std::dynamic_pointer_cast<CGM>(stateNodes[N][j]->getFactor())->getLogMass();
			odds += mass;
		} // for
	} // for

	return odds;
} // calculateEvidence()

void extractStates(const unsigned N, 
		std::map<unsigned, emdw::RVIds>& currentStates,
		std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes
		) {
	unsigned M = stateNodes[N].size();
	
	for (unsigned i = 1; i < M; i++) {
		if (stateNodes[N][i] == 0) continue;

		// Moment match the current marginal
		rcptr<Factor> marginal = stateNodes[N][i]->marginalize(elementsOfX[currentStates[N][i]], true);
		std::vector<rcptr<Factor>> comps = std::dynamic_pointer_cast<CGM>( marginal )->getComponents();
		
		// Get the mean and mass
		for (unsigned j = 0; j < comps.size(); j++) {
			double mass = std::dynamic_pointer_cast<GC>(comps[j])->getLogMass();
			ColVector<double> mean =  std::dynamic_pointer_cast<GC>(comps[j])->getMean();

			std::cout << N+1 << "," << i << "," << j << "," << mean[0] << "," << mean[2] << "," << mean[4] 
				<< "," << mass << std::endl;
		} // for
	} // for
} // extractStates()

unsigned factorial (const unsigned N) {
	return (N == 1 || N == 0) ? 1 : factorial(N-1)*N;
} // factorial()
