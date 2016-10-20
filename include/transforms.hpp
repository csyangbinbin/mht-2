/*************************************************************************
 *  Compilation: ./run_main.sh
 *  Execution: ./run_main.sh
 *  Dependencies: None
 *
 * Header file for non-linear transformations. Includes motion and
 * observation models to be used by the GaussConical class in constructing
 * joint Gaussian distributions.
 *************************************************************************/
#ifndef TRANSFORMS_HPP
#define TRANSFORMS_HPP

/**
 * This class implements the non-linear motion model function 
 * used in predicting a target's state at the next time step 
 * within a UKF. 
 *
 * @author SCJ Robertson
 * @since 11/10/16
 */
class MotionModel : public V2VTransform{
	public:
		/**
		 * Constructor for MotionModel.
		 *
		 * @param time_step The discrete time step over which the target is advanced.
		 */
		MotionModel(const double time_step) : delta_t_(time_step) {}

		/**
		 * Implements the motion model, \f$ \pmb{g} (\pmb{x}_{t-1}) \f$, in the prediction step
		 * \f$ \pmb{x}_{t} = \pmb{g} ( \pmb{x}_{t-1} ) + \pmb{\epsilon}_{t} \f$.
		 * This is used to create the joint Markovian distribution over \f$ \pmb{x}_{t-1} \f$ 
		 * and \f$ \pmb{x}_{t} \f$.
		 *
		 * @param x The previous state vector \f$ \pmb{x}_{t-1} \f$ to be transformed to \f$ \pmb{x}_{t} \f$.
		 * @return The predicted state vector \f$ \pmb{x}_{t} \f$.
		 */
		std::vector< ColVector<double> > operator()(const ColVector<double>& x) const;
	
	private:
		double delta_t_;
}; // MotionModel

/**
 * This class implements the non-linear sensor model function 
 * used in creating the observation likelihood within a UKF.
 *
 * @author SCJ Robertson
 * @since 11/10/16
 */
class SensorModel : public V2VTransform{
	public:
		/**
		 * Default constructor.
		 *
		 * @param sensor_location The position of the sensor in cartesian coordinates.
		 */
		SensorModel(const ColVector<double> &sensor_location) {
			ASSERT(sensor_location.size() == 3, "The sensor location must be defined in 3D space.");
			sensor_position_ = sensor_location;
		}
	
		/**
		 * Implements the sensor model, \f$ \pmb{h} (\pmb{x}_{t}) \f$, in the prediction step
		 * \f$ \pmb{z}_{t} = \pmb{h} ( \pmb{x}_{t} ) + \pmb{\delta}_{t} \f$.
		 * This is used to create the joint likelihood distribution over \f$ \pmb{x}_{t} \f$ 
		 * and \f$ \pmb{z}_{t} \f$.
		 *
		 * @param x The current state vector \f$ \pmb{x}_{t} \f$ to be transformed to \f$ \pmb{z}_{t} \f$.
		 * @return The measurement vector \f$ \pmb{z}_{t} \f$.
		 */
		std::vector< ColVector<double> > operator()(const ColVector<double>& x) const;

	private:
		ColVector<double> sensor_position_;
}; // SensorModel

#endif // TRANSFORMS_HPP
