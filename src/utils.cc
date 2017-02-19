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
	for (unsigned i = 0; i < L; i++) {
		globalVariables.push_back(N + i);
		map[M].push_back(N + i);
	}

	return localVariables[M];
} // addVariables()

void printCGM (const rcptr<Factor>& factor) {
	rcptr<CanonicalGaussianMixture> cgm = std::dynamic_pointer_cast<CanonicalGaussianMixture>(factor);
	std::vector<rcptr<Factor>> components = cgm->getComponents();

	std::cout << "Number of components: " << cgm->getNumberOfComponents() << "\n\n" << std::endl;
	for (rcptr<Factor> c : components) std::cout << *c << std::endl;
}
