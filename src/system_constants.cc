/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies:
 *
 * Runs intialisation code for the system constants.
 *************************************************************************/
#include "system_constants.hpp"

// Discrete time step
const double mht::kTimeStep = 0.04;

// Sensor location information
const unsigned short mht::kNumSensors = 6;
std::vector<ColVector<double>> mht::kSensorLocation; 

// Dimension of state and observation vectors
const unsigned short mht::kStateSpaceDim = 6;
const unsigned short mht::kMeasSpaceDim = 2;

// Motion model
rcptr<V2VTransform> mht::kMotionModel;
Matrix<double> mht::kRCovMat;

// Measurement model
const double mht::kC = 3e8;
const double mht::kFc = 10.525e9;
const double mht::kTp = 200e-6; 
const double mht::kBw = 47e6; 

std::vector<rcptr<V2VTransform>> mht::kMeasurementModel;
std::vector<Matrix<double>> mht::kQCovMat;

// Mahanolobis thresholding distance
const double mht::kValidationThreshold = 4;

// Smoothing paramaters
const unsigned mht::kNumberOfBackSteps = 2;

// Clutter distribution
std::vector<ColVector<double>> mht::kClutterMean;
std::vector<Matrix<double>> mht::kClutterCov;

// Gaussian mixture pruning parameters
const unsigned mht::kMaxComponents = 3;
const double mht::kThreshold = -500;
const double mht::kMergeDistance = 9;

// Launch locations
std::vector<ColVector<double>> mht::kLaunchStateMean;
std::vector<Matrix<double>> mht::kLaunchStateCov;

// Generic Means
ColVector<double> mht::kGenericMean;
Matrix<double> mht::kGenericCov;
std::vector<double> mht::kGenericWeight;

// Maximum number of targets
const unsigned mht::maxNumberOfTargets = mht::kNumSensors + 3;

// Time Off-set
const unsigned mht::timeOffSet = 0;

bool init = initialiseVariables();

// Variable management
emdw::RVIds variables;
emdw::RVIds vecX;
emdw::RVIds vecZ;
std::map<unsigned, emdw::RVIds> currentStates;
std::map<unsigned, emdw::RVIds> currentMeasurements;
std::map<unsigned, emdw::RVIds> elementsOfX;
std::map<unsigned, emdw::RVIds> elementsOfZ;
std::map<unsigned, emdw::RVIds> presentAt;
emdw::RVIds virtualMeasurementVars;

// Measurement management
rcptr<MeasurementManager> measurementManager;
unsigned kNumberOfTimeSteps;

// GraphBuilder
rcptr<GraphBuilder> graphBuilder;

// Graph representation
std::map<unsigned, std::vector<rcptr<Node>>> stateNodes;
std::map<unsigned, std::vector<rcptr<Node>>> measurementNodes;
std::vector<rcptr<Factor>> predMarginals;
std::map<unsigned, std::vector<rcptr<Factor>>> predMeasurements;
std::map<unsigned, std::vector<rcptr<Factor>>> validationRegion;

bool initialiseVariables() {
	mht::kTimeStep;

	// Sensor locations
	mht::kNumSensors;
	mht::kSensorLocation = initialiseSensorLocations();

	// Sensor location
	mht::kRCovMat = initialiseRCovMat();
	mht::kMotionModel = initialiseMotionModel();

	// Measurement models
	mht::kQCovMat = initialiseQCovMat();
	mht::kMeasurementModel = initialiseMeasurementModels();

	// Clutter cov
	mht::kClutterMean = initialiseClutterMean();
	mht::kClutterCov = initialiseClutterCovMat();
	
	// Launch covs
	mht::kLaunchStateMean = initialiseLaunchStateMean();
	mht::kLaunchStateCov = initialiseLaunchStateCov();

	// Generic covs
	mht::kGenericWeight = initialiseGenericWeights();
	mht::kGenericMean = initialiseGenericMean();
	mht::kGenericCov = initialiseGenericCov();

	return true;
} // initialiseVariables();

Matrix<double> initialiseRCovMat () {
	Matrix<double> RCov = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);

	RCov(0, 0) = 1; RCov(1, 1) = 9;
	RCov(2, 2) = 1; RCov(3, 3) = 9;
	RCov(4, 4) = 1; RCov(5, 5) = 9;

	return RCov;
} // initialiseRCovMat()

std::vector<Matrix<double>> initialiseQCovMat () {
	std::vector<Matrix<double>> QCov(mht::kNumSensors);
	
	QCov[0] = gLinear::zeros<double>(mht::kMeasSpaceDim, mht::kMeasSpaceDim);
	(QCov[0])(0, 0) = 9; (QCov[0])(1, 1) = 4;

	QCov[1] = gLinear::zeros<double>(mht::kMeasSpaceDim, mht::kMeasSpaceDim);
	(QCov[1])(0, 0) = 9; (QCov[1])(1, 1) = 4;

	QCov[2] = gLinear::zeros<double>(mht::kMeasSpaceDim, mht::kMeasSpaceDim);
	(QCov[2])(0, 0) = 9; (QCov[2])(1, 1) = 4;

	QCov[3] = gLinear::zeros<double>(mht::kMeasSpaceDim, mht::kMeasSpaceDim);
	(QCov[3])(0, 0) = 9; (QCov[3])(1, 1) = 4;

	QCov[4] = gLinear::zeros<double>(mht::kMeasSpaceDim, mht::kMeasSpaceDim);
	(QCov[4])(0, 0) = 9; (QCov[4])(1, 1) = 4;

	QCov[5] = gLinear::zeros<double>(mht::kMeasSpaceDim, mht::kMeasSpaceDim);
	(QCov[5])(0, 0) = 9; (QCov[5])(1, 1) = 4;

	return QCov;
} // initialiseQCovMat()

std::vector<ColVector<double>> initialiseClutterMean() {
	std::vector<ColVector<double>> clutterMean(mht::kNumSensors);

	// Sensor 0
	clutterMean[0] = ColVector<double>(mht::kStateSpaceDim); clutterMean[0].assignToAll(0.0);
	clutterMean[0][0] = -90; clutterMean[0][1] = 94.434;
	clutterMean[0][2] = -10; clutterMean[0][3] = -32.88;
	clutterMean[0][4] = 20; clutterMean[0][5] = 1.03;

	// Sensor 1
	clutterMean[1] = ColVector<double>(mht::kStateSpaceDim); clutterMean[1].assignToAll(0.0);
	clutterMean[1][0] = -90; clutterMean[1][1] = 98;
	clutterMean[1][2] = -10; clutterMean[1][3] = 20;
	clutterMean[1][4] = 20; clutterMean[1][5] = 0.2;

	// Sensor 2
	clutterMean[2] = ColVector<double>(mht::kStateSpaceDim); clutterMean[2].assignToAll(0.0);
	clutterMean[2][0] = -90; clutterMean[2][1] = -61;
	clutterMean[2][2] = -10; clutterMean[2][3] = -78;
	clutterMean[2][4] = 20; clutterMean[2][5] = -0.3359;

	// Sensor 3
	clutterMean[3] = ColVector<double>(mht::kStateSpaceDim); clutterMean[3].assignToAll(0.0);
	clutterMean[3][0] = -90; clutterMean[3][1] = -94;
	clutterMean[3][2] = -10; clutterMean[3][3] = 32;
	clutterMean[3][4] = 20; clutterMean[3][5] = 1.03;

	// Sensor 4
	clutterMean[4] = ColVector<double>(mht::kStateSpaceDim); clutterMean[4].assignToAll(0.0);
	clutterMean[4][0] = -90; clutterMean[4][1] = -88;
	clutterMean[4][2] = -10; clutterMean[4][3] = -45;
	clutterMean[4][4] = 20; clutterMean[4][5] = -4;

	// Sensor 5
	clutterMean[5] = ColVector<double>(mht::kStateSpaceDim); clutterMean[5].assignToAll(0.0);
	clutterMean[5][0] = -90; clutterMean[5][1] = 94;
	clutterMean[5][2] = -10; clutterMean[5][3] = 32;
	clutterMean[5][4] = 20; clutterMean[5][5] = 4;

	return clutterMean;
} // initialiseClutterMean()

std::vector<Matrix<double>> initialiseClutterCovMat () {
	std::vector<Matrix<double>> clutterCov(mht::kNumSensors);
	
	// Sensor 0
	clutterCov[0] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) (clutterCov[0])(i, i) = 625;

	// Sensor 1
	clutterCov[1] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) (clutterCov[1])(i, i) = 625;

	// Sensor 2
	clutterCov[2] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) (clutterCov[2])(i, i) = 625;

	// Sensor 3
	clutterCov[3] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) (clutterCov[3])(i, i) = 625;

	// Sensor 4
	clutterCov[4] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) (clutterCov[4])(i, i) = 625;

	// Sensor 5
	clutterCov[5] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) (clutterCov[5])(i, i) = 625;

	return clutterCov;
} // initialiseClutterCovMat()

std::vector<ColVector<double>> initialiseSensorLocations() {
	std::vector<ColVector<double>> locations(mht::kNumSensors);

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
	locations[3][0] = -185.7750; locations[3][1] = 43.614;
	locations[3][2] = -0.8730;

	// Sensor 5
	locations[4] = ColVector<double>(3); locations[4] *= 0;
	locations[4][0] = -213.8330; locations[4][1] = -47.7640;
	locations[4][2] = -8.6430;

	// Sensor 5
	locations[5] = ColVector<double>(3); locations[5] *= 0;
	locations[5][0] = 28.9590; locations[5][1] = 77.7540;
	locations[5][2] = 3.3760;

	return locations;
} // initialiseSensorLocations()

rcptr<V2VTransform> initialiseMotionModel() {
	return uniqptr<V2VTransform>(new MotionModel(mht::kTimeStep));
} // initialiseMotionModel()

std::vector<rcptr<V2VTransform>> initialiseMeasurementModels() {
	std::vector<rcptr<V2VTransform>> models(6);
	std::vector<ColVector<double>> locations = initialiseSensorLocations();

	for (unsigned i = 0; i < 6; i++) {
		models[i] = uniqptr<V2VTransform>(new SensorModel(locations[i], mht::kC, mht::kFc, mht::kTp, mht::kBw));
	}

	return models;
} // initialiseMeasurementModels()

std::vector<ColVector<double>> initialiseLaunchStateMean() {
	std::vector<ColVector<double>> launchState(3);

	// Launch state 1
	launchState[0] = ColVector<double>(mht::kStateSpaceDim); launchState[0] *= 0;
	launchState[0][0] = -12.331; launchState[0][1] = -10;
	launchState[0][2] = 13.851; launchState[0][3] = -10;
	launchState[0][4] = -5.139; launchState[0][5] = 20;

	// Launch state 2
	launchState[1] = ColVector<double>(mht::kStateSpaceDim); launchState[1] *= 0;
	launchState[1][0] = -12.997; launchState[1][1] = -10;
	launchState[1][2] = 17.325; launchState[1][3] = -10;
	launchState[1][4] = -5.151; launchState[1][5] = 20;

	// Launch state 3
	launchState[2] = ColVector<double>(mht::kStateSpaceDim); launchState[2] *= 0;
	launchState[2][0] = -13.65; launchState[2][1] = -10;
	launchState[2][2] = 20.784; launchState[2][3] = -10;
	launchState[2][4] = -5.186; launchState[2][5] = 20;

	return launchState;
} // initialiseLaunchStateMean()

std::vector<Matrix<double>> initialiseLaunchStateCov() {
	std::vector<Matrix<double>> launchCov(3);

	launchCov[0] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	launchCov[1] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	launchCov[2] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);

	for (unsigned i = 0; i < 3; i++) {
		launchCov[i](0, 0) = 1; launchCov[i](1, 1) = 9;
		launchCov[i](2, 2) = 1; launchCov[i](3, 3) = 9;
		launchCov[i](4, 4) = 1; launchCov[i](5, 5) = 9;
	}

	return launchCov;
} // initialiseLaunchStateCov()


ColVector<double> initialiseGenericMean() {
	ColVector<double> genericMean(mht::kStateSpaceDim); genericMean *= 0;
	
	genericMean[0] = -13.797; genericMean[1] = -30;
	genericMean[2] = 16.525; genericMean[3] = -15;
	genericMean[4] = -3.5824; genericMean[5] = 15;

	return genericMean;
} // initialiseGenricMean()

Matrix<double> initialiseGenericCov() {
	Matrix<double> genericCov = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	
	genericCov(0, 0) = 4; genericCov(1, 1) = 100;
	genericCov(2, 2) = 4; genericCov(3, 3) = 100;
	genericCov(4, 4) = 4; genericCov(5, 5) = 100;

	return genericCov;
} // initialiseGenericMean()

std::vector<double> initialiseGenericWeights() {
	std::vector<double> weights = {1.0/3, 1.0/3, 1.0/3};
	return weights;
} // initialiseGenericWeights()
