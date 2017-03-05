/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file declaring function prototypes for the alogrithmic
 * steps used in main.cc
 *************************************************************************/
#ifndef ALGORITHMICSTEPS_HPP
#define ALGORITHMICSTEPS_HPP

#include <iostream>
#include <map>
#include "emdw.hpp"
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "linear_gaussian.hpp"
#include "transforms.hpp"
#include "utils.hpp"
#include "system_constants.hpp"

/**
 * @brief Predict the current state of the x variables
 *
 * @param N The current time index.
 */
void predictStates(const unsigned N);

/**
 * @brief Forms hypotheses and creates current measurement distributions.
 *
 * @param N The current time index.
 */
void createMeasurementDistributions(const unsigned N);

/**
 * @brief Performs measurement update on exisitng targets.
 *
 * @param N The current time step.
 */
void measurementUpdate(const unsigned N);

#endif // ALGORITHMICSTEPS_HPP
