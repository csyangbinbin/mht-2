/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for the measurement manager class. A wrapper type which dolls
 * out sets of measurements for each sensor.
 *************************************************************************/
#ifndef MEASUREMENTMANAGER_HPP
#define MEASUREMENTMANAGER_HPP

#include <map>
#include "genvec.hpp"
#include "genmat.hpp"
#include "emdw.hpp"

/**
 *
 * @author SCJ Robertson
 * @since 16/02/17
 */
class MeasurementManager{

	public:
		/**
		 *
		 */
		MeasurementManager(const std::string& fileNameFormat, const unsigned N);

		/**
		 *
		 */
		~MeasurementManager();

	public:
		/**
		 *
		 */
		std::vector<ColVector<double>> getSensorPoints(const unsigned i, const unsigned j);

		/**
		 *
		 */
		std::vector<emdw::RVVals> getSensorMeasurements(const unsigned i, const unsigned j);

		/**
		 *
		 */
		unsigned getNumberOfTimeSteps() const;
	
	private:
		/**
		 *
		 */
		std::map<unsigned, std::map<unsigned, std::vector< ColVector<double> >>> readMeasurements(const std::string& fileNameFormat, const unsigned N);

		/**
		 *
		 */
		std::map<unsigned, std::vector< ColVector<double> >> readMeasurementFile(const std::string& fileName) const;

	private:
		std::map<unsigned, std::map<unsigned, std::vector< ColVector<double> >>> sensor_; // Fun!
		unsigned N_;
		unsigned M_;
};

#endif // MEASUREMENTMANAGER
