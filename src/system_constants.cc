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
std::vector<ColVector<double>> mht::kSensorLocation; //= initialiseSensorLocations();

// Dimension of state and observation vectors
const unsigned short mht::kStateSpaceDim = 6;
const unsigned short mht::kMeasSpaceDim = 2;

// Motion model
rcptr<V2VTransform> mht::kMotionModel = uniqptr<V2VTransform>(new MotionModel(mht::kTimeStep));
Matrix<double> mht::kRCovMat;

// Measurement model
const double mht::kC = 3e8;
const double mht::kFc = 10.525e9;
const double mht::kTp = 200e-6; 
const double mht::kBw = 47e6; 

std::vector<rcptr<V2VTransform>> mht::kMeasurementModel;
Matrix<double> mht::kQCovMat;

// Gaussian mixture pruning parameters
const unsigned mht::kMaxComponents = 100;
const double mht::kThreshold = 0.01;
const double mht::kMergeDistance = 5;

// Launch locations
std::vector<ColVector<double>> mht::kLaunchStateMean;
std::vector<Matrix<double>> mht::kLaunchStateCov;

ColVector<double> mht::kGenericMean;
Matrix<double> mht::kGenericCov;
std::vector<double> mht::kGenericWeights;

bool init = initialiseVariables();

// Variable management
emdw::RVIds variables;
emdw::RVIds vecX;
emdw::RVIds vecZ;
std::map<unsigned, emdw::RVIds> currentStates;
std::map<unsigned, emdw::RVIds> elementsOfX;
std::map<unsigned, emdw::RVIds> elementsOfZ;
std::map<unsigned, emdw::RVIds> presentAt;
emdw::RVIds virtualMeasurementVars;

// Measurement management
rcptr<MeasurementManager> manager;
unsigned kNumberOfTimeSteps;

// Graph representation
std::map<unsigned, std::vector<rcptr<Node>>> stateNodes;
std::map<unsigned, std::vector<rcptr<Node>>> measurementNodes; 
std::map<unsigned, std::vector<rcptr<Factor>>> predMeasurements;
std::map<unsigned, std::vector<rcptr<Factor>>> validationRegion;

bool initialiseVariables() {
	initialiseSensorLocations(mht::kSensorLocation);

	// Sensor location
	initialiseRCovMat(mht::kRCovMat);
	
	// Measurement models
	initialiseQCovMat(mht::kQCovMat);
	initialiseMeasurementModels(mht::kMeasurementModel, mht::kSensorLocation);
	
	// Launch covs
	initialiseLaunchStateMean(mht::kLaunchStateMean);
	initialiseLaunchStateCov(mht::kLaunchStateCov);

	// Generic covs
	initialiseGenericMean(mht::kGenericMean);
	initialiseGenericCov(mht::kGenericCov);
	initialiseGenericWeights(mht::kGenericWeights);

	return true;
} // initialiseVariables();


void initialiseRCovMat (Matrix<double>& RCov) {
	RCov = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	
	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) {
		if (i % 2 == 0) RCov(i, i) = 1;
		else RCov(i, i) = 5;
	}
} // initialiseRCovMat()

void initialiseQCovMat (Matrix<double>& QCov) {
	QCov = gLinear::zeros<double>(mht::kMeasSpaceDim, mht::kMeasSpaceDim);
	QCov(0, 0) = 9; QCov(1, 1) = 1.5;
} // initialiseQCovMat()

void initialiseSensorLocations(std::vector<ColVector<double>>& locations) {
	locations.resize(mht::kNumSensors);

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
} // initialiseSensorLocations()

void initialiseMeasurementModels(std::vector<rcptr<V2VTransform>>& models, std::vector<ColVector<double>>& locations) {
	models.resize(mht::kNumSensors);

	for (unsigned i = 0; i < 6; i++) {
		models[i] = uniqptr<V2VTransform>(new SensorModel(locations[i], 3e8, 10.525e9, 200e-6, 47e6));
	}
} // initialiseMeasurementModels()


void initialiseLaunchStateMean(std::vector<ColVector<double>>& launchState) {
	launchState.resize(3);

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
} // initialiseLaunchStateMean()

void initialiseLaunchStateCov(std::vector<Matrix<double>>& launchCov) {
	launchCov.resize(3);

	launchCov[0] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	launchCov[1] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);
	launchCov[2] = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);

	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) {
		if (i % 2 == 0) {
			(launchCov[0])(i, i) = 0.5;
			(launchCov[1])(i, i) = 0.5;
			(launchCov[2])(i, i) = 0.5;
		} else {
			(launchCov[0])(i, i) = 5;
			(launchCov[1])(i, i) = 5;
			(launchCov[2])(i, i) = 5;
		}
	}
} // initialiseLaunchStateCov()


void initialiseGenericMean(ColVector<double>& genericMean) {
	genericMean.resize(mht::kStateSpaceDim); genericMean *= 0;
	
	genericMean[0] = -12.331; genericMean[1] = -10;
	genericMean[2] = 13.851; genericMean[3] = -10;
	genericMean[4] = -5.139; genericMean[5] = 20;
} // initialiseGenricMean()

void initialiseGenericCov(Matrix<double>& genericCov) {
	genericCov = gLinear::zeros<double>(mht::kStateSpaceDim, mht::kStateSpaceDim);

	for (unsigned i = 0; i < mht::kStateSpaceDim; i++) {
		if (i % 2 == 0) (genericCov)(i, i) = 0.5;
		else (genericCov)(i, i) = 5;
	}
} // initialiseGenericMean()

void initialiseGenericWeights(std::vector<double>& weights) {
	weights.clear(); weights = {1.0};
} // initialiseGenericWeights()
