/*************************************************************************
 *  Compilation:  
 *  Execution:    
 *  Dependencies:
 *
 * Main app, runs small examples for now.
 *************************************************************************/

// patrec headers
#include "error.hpp"
#include "testing.hpp"
#include "genmat.hpp"
#include "vecset.hpp"
#include "mean_cov.hpp"

// emdw headers
#include "emdw.hpp"
#include "anytype.hpp"
#include "clustergraph.hpp"
#include "gausscanonical.hpp"
#include "oddsandsods.hpp"
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
 * Creates a random with each element drawn independently
 * from a random a Gaussian distribution.
 *
 * @param size_x The dimension of the column vector.
 * @param generator The rquired random number generator.
 * @param mean The mean of the Gaussian distribution.
 * @param var The variance of the Gaussian distribution.
 * @return A Gaussian random vector.
 */
ColVector<double> NormVec(int size_x, default_random_engine generator, double mean, double var) {
	
	ColVector<double> randvec(size_x);
	normal_distribution<double> dist(mean, sqrt(var));

	for (int i = 0; i < size_x; i++) randvec[i] = dist(generator);
	
	return randvec;
}

/**
 * Reads the mean from a factor, rather than acessing a member function. For some
 * reason. Seems to specifically apply to a Cannonical representation.
 *
 * @param factor The factor, supposedly in cannonical form.
 * @return The mean of the factor.
 */
ColVector<double> ReadMean(rcptr<Factor> factor) {
	ColVector<double> mean;
	Matrix<double> cov;

	emdw::RVIds vars;
	Matrix<double> K_mat;
	ColVector<double> h_vec;

	stringstream file;
	file << *factor;
	string dummy;

	file >> dummy;
	CHECK(file >> vars, IOError, "Error reading variable's RowVector");

	file >> dummy;
	CHECK(file >> K_mat, IOError, "Error reading K.");

	file >> dummy;
	CHECK(file >> h_vec, IOError, "Error reading h.");

	int fail;
	cov = inv(K_mat, fail);

	mean = cov*h_vec;
	return mean;
}

/**
 * Find the maximum difference between elements in two vectors of equal length.
 *
 * @param vector_1 The first vector.
 * @param vector_2 The second vector, must be equal in length to vector_1
 * @return The maximum difference between elements.
 */
double MaxDiff(ColVector<double> vector_1, ColVector<double> vector_2) {
	double max_dif = 0;
	double dif;

	for (int i = 0; i < vector_1.size(); i++) {
		dif = fabs(vector_1[i] - vector_2[i])/min(vector_1[i], vector_2[i]);
		if (dif > max_dif) max_dif = dif;
	}

	return max_dif;
}

/**
 * Creates a map between the observed variables and their observed states. 
 *
 * @param observed_vars The variables to be fixed in their observed state.
 * @param observed_vals The variables' repsective observed states.
 * @return A map of the variables and their given observed states.
 */
std::map<emdw::RVIdType, AnyType> MapObserved(const emdw::RVIds& observed_vars, const emdw::RVVals& observed_vals) { 
	
	std::map<emdw::RVIdType, AnyType> observed;

	for (unsigned i = 0; i < observed_vars.size(); i++) {
		observed[observed_vars[i]] = observed_vals[i];
	}

	return observed;
}

/**
 * Main app, runs small examples for now.
 *
 * @author SCJ Robertson
 * @since 0.9.1
 */
int main(int, char *argv[]) {
	typedef GaussCanonical GCT1;
	
	const int no_obs = 5;
	double seed = 100;
	default_random_engine generator(seed);

	double x_var = 100;
	double r_var = 0.1*1.5;
	double q_var = 0.1*2;

	Matrix<double> A_mat(3, 3);
	A_mat(0, 0) = 3; A_mat(0, 1) = 2; A_mat(0, 2) = 1;
	A_mat(1, 0) = 2; A_mat(1, 1) = 4; A_mat(1, 2) = 1;
	A_mat(2, 0) = 1; A_mat(2, 1) = 1; A_mat(2, 2) = 3;

	Matrix<double> C_mat(3, 3);
	C_mat(0, 0) = 1; C_mat(0, 1) = 3; C_mat(0, 2) = 3;
	C_mat(1, 0) = 3; C_mat(1, 1) = 5; C_mat(1, 2) = 1;
	C_mat(2, 0) = 3; C_mat(2, 1) = 1; C_mat(2, 2) = 3;

	Matrix<double> R_mat = gLinear::zeros<double>(3, 3);
	R_mat(0, 0) = r_var; R_mat(1, 1) = r_var; R_mat(2, 2) = r_var;

	Matrix<double> Q_mat = gLinear::zeros<double>(3, 3);
	Q_mat(0, 0) = q_var; Q_mat(1, 1) = q_var; Q_mat(2, 2) = q_var;
	
	ColVector<double> mu_0(3);
	mu_0[0] = 1; mu_0[1] = 2; mu_0[2] = 3;

	vector<ColVector<double>> x_vals;
	x_vals.push_back(mu_0 + NormVec(3, generator, 0, x_var));

	vector<ColVector<double>> x_naive;
	x_naive.push_back(mu_0);

	for (int i = 0; i < (no_obs -1); i++) {
		x_vals.push_back(A_mat*x_vals[i] + NormVec(3, generator, 0, r_var));
		x_naive.push_back(A_mat*x_naive[i]);
	}
	
	vector<ColVector<double>> y_vals;
	for (unsigned i = 0; i < x_vals.size(); i++) {
		y_vals.push_back(C_mat*x_vals[i] + NormVec(3, generator, 0, q_var));
	}

	RVIds vars = {0, 1, 2};
	ColVector<double> mean = mu_0;
	Matrix<double> cov = gLinear::zeros<double>(3, 3);
	cov(0, 0) = x_var; cov(1, 1) = x_var; cov(2, 2) = x_var;
	
	int fail;
	Matrix<double> K_mat = inv(cov, fail);
	ColVector<double> h_vec = K_mat*mean;

	rcptr<Factor> gauss_x0 = uniqptr<GCT1> ( new GCT1(vars, K_mat, h_vec, true) );
	vector<rcptr<Factor>> factors;
	factors.push_back(gauss_x0);

	RVIds prev_vars = {0, 1, 2};
	for (int i = 0; i < (no_obs-1); i++) {
		unsigned v = (i+1)*10;
		vars = {v, v+1, v+2};

		factors.push_back( uniqptr<GCT1>( new GCT1 (factors[i]->marginalize(prev_vars)->copy(), A_mat, vars, R_mat) ) );

		prev_vars = {v, v+1, v+2};
	}

	int size = factors.size();
	prev_vars = {0, 1, 2};
	for (int i = 0; i < size; i++) {
		unsigned v = i*10;
		vars = {v+3, v+4, v+5};

		factors.push_back( uniqptr<GCT1>( new GCT1 (factors[i]->marginalize(prev_vars)->copy(), C_mat, vars, Q_mat) ) );
		
		v = (i+1)*10;
		prev_vars = {v, v+1, v+2};
	}

	RVVals values;
	
	map<Idx2, rcptr<Factor>> messages;
	MessageQueue msg_queue;

	values.clear();
	vars.clear();
	for (unsigned i = no_obs; i < factors.size(); i++) {
		rcptr<Factor> temp_ptr;

		for(int j = 0; j < 3; j++) values.push_back(y_vals[i-no_obs][j]);

		unsigned v = (i-no_obs)*10;
		vars.push_back(v+3);
		vars.push_back(v+4);
		vars.push_back(v+5);
	}

	ClusterGraph hcg_loopy(factors, MapObserved(vars, values));
	hcg_loopy.exportToGraphViz("kalman");

	loopyBP_CG(hcg_loopy, messages, msg_queue, 0.0);

	typedef map<Idx2, rcptr<Factor>>::iterator it_type;
	it_type iterator = messages.begin();

	rcptr<Factor> temp_ptr1 = iterator->second;
	iterator++;
	rcptr<Factor> temp_ptr2 = iterator->second;

	ColVector<double> temp_mean;
	double tolerance = 0.05;

	for(it_type iterator = messages.begin(); iterator != messages.end(); iterator++) {
		
		rcptr<Factor> temp_ptr1 = iterator->second;
		iterator++;
		rcptr<Factor> temp_ptr2 = iterator->second;

		cout << *(temp_ptr1->absorb(temp_ptr2)) << endl;
		temp_mean = ReadMean((temp_ptr1->absorb(temp_ptr2)));
	}
	
	
	/*
	for(unsigned i = no_obs; i < factors.size(); i++) {
		rcptr<Factor> temp_ptr;
		values.clear();

		for (int j = 0; j < 3; j++) values.push_back(y_vals[i-no_obs][j]);

		unsigned v = (i-no_obs)*10;
		vars = {v+3, v+4, v+5};

		temp_ptr = factors[i]->observeAndReduce(vars, values);
		factors[i] = rcptr<Factor>(temp_ptr->copy());		
	}


	factors[1]->inplaceAbsorb(factors[no_obs]);
	for (unsigned i = no_obs; i < (factors.size()-1); i++) {
		factors[i-no_obs]->inplaceAbsorb(factors[i+1]);
	}

	vector<rcptr<Factor>> messages;
	vars = {10, 11, 12};
	messages.push_back(rcptr<Factor> (factors[1]->marginalize(vars)->copy()));

	for (unsigned i = 2; i < (no_obs -1); i++) {
		unsigned v = i*10;

		vars = {v, v+1, v+2};
		rcptr<Factor> temp;
		temp = messages[i-2]->absorb(factors[i]);
		temp = temp->marginalize(vars);

		messages.push_back(rcptr<Factor>(temp->copy()));
	}

	unsigned v = (no_obs-2)*10;
	vars = {v, v+1, v+2};
	messages.push_back(rcptr<Factor>(factors[no_obs-1]->marginalize(vars)->copy()));

	for (unsigned i = (no_obs-3), j = (no_obs-2); i > 0; i--, j++) {
		unsigned v = i*10;
		vars = {v, v+1, v+2};

		rcptr<Factor> temp;
		temp = messages[j]->absorb(factors[i+1]);
		temp = temp->marginalize(vars);

		messages.push_back( rcptr<Factor>(temp->copy()) );
	}

	double tolerance = 0.05;
	int index = messages.size()-1;
	ColVector<double> temp_mean = ReadMean(factors[1]->absorb(messages[index])->marginalize(RVIds{0, 1, 2}));

	if (MaxDiff(temp_mean, x_vals[0]) > tolerance) {
		LOGERROR(pLog, "Testing", "Expected: " << x_vals[0] << endl << "Got: " << temp_mean);
		LOGERROR(pLog, "Testing", "GCF Failed");
		THROW(TestFailed, "");
	}
	*/

}
