/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file declaring system wide constants.
 *************************************************************************/
#ifndef SYSTEMCONSTANTS_HPP
#define SYSTEMCONSTANTS_HPP

#include <iostream>
#include <map>
#include "emdw.hpp"
#include "discretetable.hpp"
#include "gausscanonical.hpp"
#include "canonical_gaussian_mixture.hpp"
#include "linear_gaussian.hpp"
#include "transforms.hpp"

// Function prototypes
std::vector<ColVector<double>> initialiseSensorLocations();
std::vector<rcptr<V2VTransform>> initialiseMeasurementModels();
std::vector<ColVector<double>> initialiseLaunchStateMean();
std::vector<Matrix<double>> initialiseLaunchStateCov();

// Discrete time step
extern const double kTimeStep;

// Sensor location information
extern const unsigned short kNumSensors;
extern std::vector<ColVector<double>> kSensorLocation;

// Dimension of state and observation vectors
extern const unsigned short N;
extern const unsigned short M;

// Motion model
extern rcptr<V2VTransform> kMotionModel;

// Measurement model
extern std::vector<rcptr<V2VTransform>> kMeasurementModel;

// Gaussian mixture pruning parameters
extern const unsigned kMaxComponents;
extern const double kThreshold;
extern const double kMergeDistance;

// Launch locations
extern std::vector<ColVector<double>> kLaunchStateMean;
extern std::vector<Matrix<double>> kLaunchStateCov;

// Variable management
extern emdw::RVIds variables;
extern emdw::RVIds vectorVars;
extern std::vector<std::map<unsigned, bool>> present;
extern std::map<unsigned, emdw::RVIds> vecElements;

#endif // SYSTEMCONSTANTS_HPP
