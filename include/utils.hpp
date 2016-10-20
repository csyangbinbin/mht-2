/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * A library of useful helper functions which don't really fit in anywhere 
 * else. 
 *************************************************************************/
#ifndef UTILS_HPP
#define UTILS_HPP

/** 
 * Almost all operator return a pointer Factor rather than the derived class,
 * making it difficult to access derived class members. This function
 * returns the mean of a GaussianCanonical factor when represented as 
 * its supertype.
 *
 * @param factor A GaussCanonical factor represented in its abstract type.
 * @return The mean of the given factor.
 */
ColVector<double> ReadMean(rcptr<Factor> factor);

/** 
 * Almost all operator return a pointer Factor rather than the derived class,
 * making it difficult to access derived class members. This function
 * returns the covariance of a GaussianCanonical factor when represented as 
 * its supertype.
 *
 * @param factor A GaussCanonical factor represented in its abstract type.
 * @return The covariance of the given factor.
 */
Matrix<double> ReadCovariance(rcptr<Factor> factor);

#endif // UTILS_HPP
