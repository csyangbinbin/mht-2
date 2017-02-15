/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs intialisation code for the system constants.
 *************************************************************************/
#include "system_constants.hpp"

// Discrete time step
const double kTimeStep = 0.04;

// Sensor location information
const unsigned short kNumSensors = 6;
std::vector<ColVector<double>> kSensorLocation = initialiseSensorLocations();

// Dimension of state and observation vectors
const unsigned short N = 6;
const unsigned short M = 2;

// Motion model
rcptr<V2VTransform> kMotionModel = uniqptr<V2VTransform>(new MotionModel(kTimeStep));

// Measurement model
std::vector<rcptr<V2VTransform>> kMeasurementModel = initialiseMeasurementModels();

// Gaussian mixture pruning parameters
const unsigned kMaxComponents = 100;
const double kThreshold = 0.01;
const double kMergeDistance = 5;

// Launch locations
std::vector<ColVector<double>> kLaunchStateMean = initialiseLaunchStateMean();
std::vector<Matrix<double>> kLaunchStateCov = initialiseLaunchStateCov();
ColVector<double> kGenericMean = initialiseGenericMean();
Matrix<double> kGenericCov = initialiseGenericCov();

// Variable management
emdw::RVIds variables;
emdw::RVIds vecX;
emdw::RVIds vecZ;
std::map<unsigned, emdw::RVIds> elementsOf;

std::vector<ColVector<double>> initialiseSensorLocations() {
	std::vector<ColVector<double>> locations(kNumSensors);

	// Sensor 1
	locations[0] = ColVector<double>(3); locations[0] *= 0;
	locations[0][0] = 9.7010; locations[0][1] = -24.4460;
	locations[0][2] = 1.2630;

	// Sensor 2
	locations[1] = ColVector<double>(3); locations[1] *= 0;
	locations[1][0] = -10.0950; locations[1][1] = 37.499;
	locations[1][2] = -5.6910;

	// Sensor 3
	locations[2] = ColVector<double>(3); locations[2] *= 0;
	locations[2][0] = -173.4660; locations[2][1] = -96.6460;
	locations[2][2] = -6.67;

	// Sensor 4
	locations[3] = ColVector<double>(3); locations[3] *= 0;
	locations[3][0] = -185.7750; locations[3][1] = -43.614;
	locations[3][2] = -0.8730;

	// Sensor 5
	locations[4] = ColVector<double>(3); locations[4] *= 0;
	locations[4][0] = -213.8330; locations[4][1] = -47.7640;
	locations[4][2] = -8.6430;

	// Sensor 5
	locations[5] = ColVector<double>(3); locations[5] *= 0;
	locations[5][0] = 28.9590; locations[5][1] = 3.3760;
	locations[5][2] = 77.7540;

	return locations;
} // initialiseSensorLocations()

std::vector<rcptr<V2VTransform>> initialiseMeasurementModels() {
	std::vector<rcptr<V2VTransform>> models(kNumSensors);
	std::vector<ColVector<double>> locations = initialiseSensorLocations();

	for (unsigned i = 0; i < kNumSensors; i++) {
		models[i] = uniqptr<V2VTransform>(new SensorModel(locations[i]));
	}
		
	return models;
} // initialiseMeasurementModels()

std::vector<ColVector<double>> initialiseLaunchStateMean() {
	std::vector<ColVector<double>> launchState(3);

	// Launch state 1
	launchState[0] = ColVector<double>(N); launchState[0] *= 0;
	launchState[0][0] = -12.331; launchState[0][1] = -10;
	launchState[0][2] = 13.851; launchState[0][3] = -10;
	launchState[0][4] = -5.139; launchState[0][5] = 20;

	// Launch state 2
	launchState[1] = ColVector<double>(N); launchState[1] *= 0;
	launchState[1][0] = -12.997; launchState[1][1] = -10;
	launchState[1][2] = 17.325; launchState[1][3] = -10;
	launchState[1][4] = -5.151; launchState[1][5] = 20;

	// Launch state 3
	launchState[2] = ColVector<double>(N); launchState[2] *= 0;
	launchState[2][0] = -13.65; launchState[2][1] = -10;
	launchState[2][2] = 20.784; launchState[2][3] = -10;
	launchState[2][4] = -5.186; launchState[2][5] = 20;

	return launchState;
} // initialiseLaunchStateMean()

std::vector<Matrix<double>> initialiseLaunchStateCov() {
	std::vector<Matrix<double>> launchState(3);

	launchState[0] = gLinear::zeros<double>(N, N);
	launchState[1] = gLinear::zeros<double>(N, N);
	launchState[2] = gLinear::zeros<double>(N, N);

	for (unsigned i = 0; i < N; i++) {
		if (i % 2 == 0) {
			(launchState[0])(i, i) = 0.5;
			(launchState[1])(i, i) = 0.5;
			(launchState[2])(i, i) = 0.5;
		} else {
			(launchState[0])(i, i) = 5;
			(launchState[1])(i, i) = 5;
			(launchState[2])(i, i) = 5;
		}
	}

	return launchState;
} // initialiseLaunchStateCov()

ColVector<double> initialiseGenericMean() {
	ColVector<double> genericMean(N); genericMean *= 0;

	genericMean[0] = -14.3310; genericMean[1] = -10.00;
	genericMean[2] = 11.8510; genericMean[3] = -10.00;
	genericMean[4] = -1.3352; genericMean[5] = 18.0380;

	return genericMean;
} // initialiseGenricMean()

Matrix<double> initialiseGenericCov() {
	Matrix<double> genericCov = gLinear::zeros<double>(N, N);

	for (unsigned i = 0; i < N; i++) {
		if (i % 2 == 0) (genericCov)(i, i) = 6.440;
		else {
			(genericCov)(i, i) = 30.00;
			(genericCov)(i - 1, i) = (genericCov)(i, i - 1) = 3.0;
		} 
	}

	return genericCov;
} // initialiseGenricMean()
