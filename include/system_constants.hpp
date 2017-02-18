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
#include "node.hpp"
#include "graph.hpp"
#include "graph_builder.hpp"
#include "measurement_manager.hpp"
#include "transforms.hpp"

// Function prototypes
Matrix<double> initialiseRCovMat();
Matrix<double> initialiseQCovMat();
std::vector<ColVector<double>> initialiseSensorLocations();
std::vector<rcptr<V2VTransform>> initialiseMeasurementModels();
std::vector<ColVector<double>> initialiseLaunchStateMean();
std::vector<Matrix<double>> initialiseLaunchStateCov();
ColVector<double> initialiseGenericMean();
Matrix<double> initialiseGenericCov();

// Typedefs
typedef DiscreteTable<unsigned> DT;
typedef GaussCanonical GC;
typedef CanonicalGaussianMixture CGM;
typedef LinearGaussian LG;

// Discrete time step
extern const double kTimeStep;

// Sensor location information
extern const unsigned short kNumSensors;
extern std::vector<ColVector<double>> kSensorLocation;

// Dimension of state and observation vectors
extern const unsigned short kStateSpaceDim;
extern const unsigned short kMeasSpaceDim;

// Motion model
extern rcptr<V2VTransform> kMotionModel;
extern Matrix<double> kRCovMat;

// Measurement model
extern std::vector<rcptr<V2VTransform>> kMeasurementModel;
extern Matrix<double> kQCovMat;

// Gaussian mixture pruning parameters
extern const unsigned kMaxComponents;
extern const double kThreshold;
extern const double kMergeDistance;

// Launch locations
extern std::vector<ColVector<double>> kLaunchStateMean;
extern std::vector<Matrix<double>> kLaunchStateCov;
extern ColVector<double> kGenericMean;
extern Matrix<double> kGenericCov;

// Variable management
extern emdw::RVIds variables; // Global variables
extern emdw::RVIds vecX;
extern emdw::RVIds vecZ;
extern std::map<unsigned, emdw::RVIds> currentStates;
extern std::map<unsigned, emdw::RVIds> elementsOf;
extern std::map<unsigned, emdw::RVIds> presentAt;

// Measurement management
extern rcptr<MeasurementManager> manager;
extern unsigned kNumberOfTimeSteps;

// Graph representation
extern std::map<unsigned, std::vector<rcptr<Node>>> stateNodes;
extern std::map<unsigned, std::vector<rcptr<Node>>> measurementNodes;
extern std::map<unsigned, std::vector<rcptr<Factor>>> factors;

#endif // SYSTEMCONSTANTS_HPP
