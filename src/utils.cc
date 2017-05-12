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

double calculateEvidence(const unsigned N, std::map<unsigned, std::vector<rcptr<Node>>>& stateNodes) {
	unsigned M = stateNodes[N].size();
	double odds = 0;

	// Determine the log-odds - excluding vacuous sponge.
	for (unsigned i = 0; i < M; i++) {
		if (stateNodes[N][i] == 0) continue;
		double mass = std::dynamic_pointer_cast<CGM>(stateNodes[N][i]->getFactor())->getLogMass();
		if (N == 13) {
			std::cout << "Target " << stateNodes[N][i]->getIdentity() 
				<< ", Mass: " << mass << " : " << exp(mass) << std::endl;
		}
		odds += mass;
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
		rcptr<Factor> marginal = stateNodes[N][i]->marginalize(elementsOfX[currentStates[N][i]], true);
		rcptr<Factor> matched = std::dynamic_pointer_cast<CGM>( marginal )->momentMatch();
	
		ColVector<double> mean =  std::dynamic_pointer_cast<GC>(matched)->getMean();
		std::cout << N << ";" << i << ";" << mean[0] << ";" << mean[2] << ";" << mean[4] << std::endl;

		if (std::isnan(mean[0])) {
			std::cout << *stateNodes[N][i] << std::endl;
		}

	} // for
} // extractStates()

unsigned factorial (const unsigned N) {
	return (N == 1 || N == 0) ? 1 : factorial(N-1)*N;
} // factorial()
