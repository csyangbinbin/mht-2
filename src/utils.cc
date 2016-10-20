/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * A library of useful helper functions which don't really fit in anywhere 
 * else. 
 *************************************************************************/

#include <iostream>
#include <string>
#include "genmat.hpp"
#include "emdw.hpp"
#include "factor.hpp"
#include "utils.hpp"

ColVector<double> ReadMean(rcptr<Factor> factor) {
	emdw::RVIds variables;
	ColVector<double> mu, h;
	Matrix<double> S_mat, K_mat;

	std::stringstream file;
	std::string temp;
	
	file << *factor;
	file >> temp;
	CHECK(file >> variables, IOError, "Error reading variables.");

	file >> temp;
	CHECK(file >> K_mat, IOError, "Error reading K_mat");

	file >> temp;
	CHECK(file >> h, IOError, "Error reading h.");

	int failure;
	S_mat = inv(K_mat, failure);
	mu = S_mat*h;
	
	return mu;
} // ReadMean()

Matrix<double> ReadCovariance(rcptr<Factor> factor) {
	emdw::RVIds variables;
	ColVector<double> mu, h;
	Matrix<double> S_mat, K_mat;

	std::stringstream file;
	std::string temp;
	
	file << *factor;
	file >> temp;
	CHECK(file >> variables, IOError, "Error reading variables.");

	file >> temp;
	CHECK(file >> K_mat, IOError, "Error reading K_mat");

	file >> temp;
	CHECK(file >> h, IOError, "Error reading h.");

	int failure;
	S_mat = inv(K_mat, failure);
	
	return S_mat;
} // ReadCovariance()

