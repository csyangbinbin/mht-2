/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Source file for the measurement manager class.
 *************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"
#include "measurement_manager.hpp"

MeasurementManager::MeasurementManager(const std::string& fileNameFormat, const unsigned N) : N_(N) {
	sensor_ = readMeasurements(fileNameFormat, N);
	M_ = sensor_[0].size();
} // Default constructor

MeasurementManager::~MeasurementManager() {  } // Default destructor


std::vector<ColVector<double>> MeasurementManager::getSensorPoints(const unsigned i, const unsigned j) const {
	return (sensor_[i])[j];
} // getSensorMeasurements()

std::vector<emdw::RVVals> MeasurementManager::getSensorMeasurements(const unsigned i, const unsigned j) const {
	std::vector<emdw::RVVals> vals;
	std::vector<ColVector<double>> points = getSensorPoints(i, j);

	for (ColVector<double> i : points) vals.push_back(emdw::RVVals{ i[0], i[1] });
	
	return vals;
} // getSensorMeasurements()

unsigned MeasurementManager::getNumberOfTimeSteps() const {
	return M_;
} // getNumberOfTimeSteps()

std::map<unsigned, std::map<unsigned, std::vector< ColVector<double> > >> MeasurementManager::readMeasurements(const std::string& fileNameFormat, const unsigned N) {
	std::map<unsigned, std::map<unsigned, std::vector< ColVector<double> >>> map;

	for (unsigned i = 0; i < N; i++) {
		std::string fileName = fileNameFormat + "_" + std::to_string(i + 1) + ".csv";
		map[i] = readMeasurementFile(fileName);
	}
	return map;
} // readMeasurements()

std::map<unsigned, std::vector< ColVector<double>  >> MeasurementManager::readMeasurementFile(const std::string& fileName) const {
	std::map<unsigned, std::vector< ColVector<double>  >> map;

	std::ifstream data(fileName);
	std::string line;

	while(std::getline(data, line)) {
		std::stringstream lineStream(line);
		std::string cell;
		std::vector<double> result;
		
		while(std::getline(lineStream, cell, ';')) result.push_back( strtod(cell.c_str(), NULL) );
		
		unsigned N = (unsigned) result[0];
		map[N].push_back( {} );

		if (std::isfinite(result[1])) 
		{
			ColVector<double> vals = ColVector<double>(2);
			vals[0] = result[1]; vals[1] = result[2];
			map[N].push_back(vals);
		} 
	}

	data.close();
	
	return map;
} // readMeasurementFile()
