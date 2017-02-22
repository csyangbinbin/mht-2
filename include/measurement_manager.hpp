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
		 * @brief Constructor
		 *
		 * Reads in measurements from a set of input files.
		 *
		 * @param fileNameFormat The location of the files and their root
		 * format.
		 *
		 * @param N The number of files in the batch.
		 */
		MeasurementManager(const std::string& fileNameFormat, const unsigned N);

		/**
		 * @brief Default destructor.
		 */
		~MeasurementManager();

	public:
		/**
		 * @brief Return a set of measurements from a sensor at a given time step.
		 *
		 * Returns a set of measurements from a sensor at a given time step. This 
		 * returns a set of gLinear::ColVectors. Missed detections are
		 * given by an empty vector.
		 *
		 * @param i The sensor number
		 *
		 * @param j The current time step indice. 
		 *
		 * @return A vector of ColVectors of Range-Doppler measurements from a particular
		 * sensor.
		 */
		std::vector<ColVector<double>> getSensorPoints(const unsigned i, const unsigned j) const;

		/**
		 * @brief Return a set of measurements from a sensor at a given time step. 
		 *
		 * Returns a set of measurements from a sensor at a given time step. This returns
		 * a set of emdw::RVVals to intorduced as evidence in a factor. Missed detections
		 * are given by an empty set.
		 *
		 * @param i The sensor number
		 *
		 * @param j The current time step indice. 
		 *
		 * @return A vector emdw::RVVals of Range-Doppler measureemtns from a particular
		 * sensor.
		 */
		std::vector<emdw::RVVals> getSensorMeasurements(const unsigned i, const unsigned j) const;

		/**
		 * @brief Returns the total number of time steps over which the data was collected.
		 *
		 * @return The total number of time steps.
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
		mutable std::map<unsigned, std::map<unsigned, std::vector< ColVector<double> >>> sensor_; // Fun!
		unsigned N_;
		unsigned M_;
};

#endif // MEASUREMENTMANAGER
