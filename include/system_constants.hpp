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
#include "conditional_gaussian.hpp"
#include "node.hpp"
#include "graph_builder.hpp"
#include "measurement_manager.hpp"
#include "transforms.hpp"

// Function prototypes
Matrix<double> initialiseRCovMat();
Matrix<double> initialiseQCovMat();
ColVector<double> initialiseClutterMean();
Matrix<double> initialiseClutterCovMat();
std::vector<ColVector<double>> initialiseSensorLocations();
rcptr<V2VTransform> initialiseMotionModel();
std::vector<rcptr<V2VTransform>> initialiseMeasurementModels();
std::vector<ColVector<double>> initialiseLaunchStateMean();
std::vector<Matrix<double>> initialiseLaunchStateCov();
ColVector<double> initialiseGenericMean();
Matrix<double> initialiseGenericCov();
std::vector<double> initialiseGenericWeights();
bool initialiseVariables();

// Typedefs
typedef unsigned short T;
typedef std::vector<T> DASS;
typedef DiscreteTable<unsigned> DT;
typedef GaussCanonical GC;
typedef CanonicalGaussianMixture CGM;
typedef ConditionalGaussian CLG;

namespace mht {
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
	extern const double kC; // Speed of light
	extern const double kFc; // Carrier frequency
	extern const double kTp; // Sweep period
	extern const double kBw; // Bandwidth

	extern std::vector<rcptr<V2VTransform>> kMeasurementModel;
	extern Matrix<double> kQCovMat;

	// Smoothing parameters
	extern const unsigned kNumberOfBackSteps;

	// Mahanalobis thresholding distance
	extern const double kValidationThreshold;

	// Clutter distribution
	extern ColVector<double> kClutterMean;
	extern Matrix<double> kClutterCov;

	// Gaussian mixture pruning parameters
	extern const unsigned kMaxComponents;
	extern const double kThreshold;
	extern const double kMergeDistance;

	// Launch locations
	extern std::vector<ColVector<double>> kLaunchStateMean;
	extern std::vector<Matrix<double>> kLaunchStateCov;

	extern ColVector<double> kGenericMean;
	extern Matrix<double> kGenericCov;
	extern std::vector<double> kGenericWeight;

	// Force initialisation
	extern bool init;
} // class

// Variable management
extern emdw::RVIds variables; // Global variables
extern emdw::RVIds vecX;
extern emdw::RVIds vecZ;
extern std::map<unsigned, emdw::RVIds> currentStates;
extern std::map<unsigned, emdw::RVIds> currentMeasurements;
extern std::map<unsigned, emdw::RVIds> elementsOfX;
extern std::map<unsigned, emdw::RVIds> elementsOfZ;
extern std::map<unsigned, emdw::RVIds> presentAt;
extern emdw::RVIds virtualMeasurementVars;

// Measurement management
extern rcptr<MeasurementManager> measurementManager;
extern unsigned kNumberOfTimeSteps;

// GraphBuilder
extern rcptr<GraphBuilder> graphBuilder;

// Graph representation
extern std::map<unsigned, std::vector<rcptr<Node>>> stateNodes;
extern std::map<unsigned, std::vector<rcptr<Node>>> measurementNodes;
extern std::vector<rcptr<Factor>> predMarginals;
extern std::map<unsigned, std::vector<rcptr<Factor>>> predMeasurements;
extern std::map<unsigned, std::vector<rcptr<Factor>>> validationRegion;

#endif // SYSTEMCONSTANTS_HPP
