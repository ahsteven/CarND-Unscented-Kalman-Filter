#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  // Predicted state
  VectorXd x_pred_;

  ///* state covariance matrix
  MatrixXd P_;
  ///*predicted state covariance matrix
  MatrixXd P_pred_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

 
  

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;



  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;



  //// state transistion matrix
  //Eigen::MatrixXd F_;

  //// process covariance matrix
  //Eigen::MatrixXd Q_;


  // measurement covariance matrix
  Eigen::MatrixXd R_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);
  void GenerateSigmaPoints();
  void AugmentedSigmaPoints();
  void SigmaPointPrediction(double dt);
  void PredictMeanAndCovariance();
  //void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out);
  //void PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out);
  //void UpdateState(VectorXd* x_out, MatrixXd* P_out);


  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);

  private:
	  // check whether the tracking toolbox was initiallized or not (first measurement)
	  bool is_initialized_;
	  // previous timestamp
	  long long previous_timestamp_;
	  // tool object used to compute Jacobian and RMSE
	  Tools tools;
	  Eigen::MatrixXd R_laser_;
	  Eigen::MatrixXd R_radar_;
	  Eigen::MatrixXd H_laser_;
	  VectorXd z_out_ ;
	  MatrixXd S_out_ ;
	  MatrixXd Xsig;

	  //create augmented mean vector
	  VectorXd x_aug;
	  //create augmented state covariance
	  MatrixXd P_aug;
	  //create sigma point matrix
	  MatrixXd Xsig_aug;
	  VectorXd weights_;
	  double lambda_aug_;// spreading parameter for augmented state
	  double lambda_;// spreading parameter
					 //create matrix for sigma points in measurement space
	  // Holds Sigma points in radar measuremente space
	  MatrixXd Zsig_R_;
	  // Holds Sigma points in lidar measuremente space
	  MatrixXd Zsig_L_;
	  //create matrix for cross correlation Tc for radar
	  MatrixXd Tc_R_;
	  //create matrix for cross correlation Tc for lidar
	  MatrixXd Tc_L_;
};

#endif /* UKF_H */
