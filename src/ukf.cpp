#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // First run
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);
  x_pred_ = VectorXd(n_x_);
  // initial covariance matrix
  //P_ = MatrixXd(n_x_, n_x_);
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_pred_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = .9;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  R_radar_ = MatrixXd(3, 3);
  R_laser_ = MatrixXd(2, 2);
  H_laser_ = MatrixXd(2, 5);
  H_laser_ << 1, 0, 0, 0, 0,
	          0, 1, 0, 0, 0;
  // sigma points initialize
  Xsig = MatrixXd(n_x_, 2 * n_x_+1);

  Zsig_R_ = MatrixXd(3, 2 * n_aug_ + 1);
  //create matrix for sigma points in measurement space
  Zsig_L_ = MatrixXd(2, 2 * n_aug_ + 1);

  //augmented mean vector
   x_aug = VectorXd(n_aug_);
  //augmented state covariance
  P_aug = MatrixXd(n_aug_, n_aug_);
  //create sigma point matrix
  Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  //create predictedsigma point matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // initialize weights for augmented mean calculation
  weights_ = VectorXd(2*n_aug_+1);
  // Cross correlation matrix radar update step
  Tc_R_ = MatrixXd(n_x_, 3);
  // Cross correlation matrix lidar update step
  Tc_L_ = MatrixXd(n_x_, 2);



}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {


	if (!is_initialized_) {



		// Initialize State
		if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
			//std::cout << "radar init" << endl;
			/**
			Convert radar from polar to cartesian coordinates and initialize state.
			*/

			float rho = meas_package.raw_measurements_(0);
			float phi = meas_package.raw_measurements_(1);
			float rhodot = meas_package.raw_measurements_(2);

			// set initial state
			//x_ << rho*cos(phi), rho*sin(phi), rhodot,phi,0;
			x_ << rho*cos(phi), rho*sin(phi), 0,0,0;
		}
		else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

			//std::cout << "laser init" << endl;
			x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
		}
		//Set spreading parameter
		lambda_ = 3 - n_x_;
		//Set augmented spreading parameter
		lambda_aug_ = 3 - n_aug_;
		// set weights
		double weight_0 = lambda_aug_/(lambda_aug_+n_aug_);
		weights_(0) = weight_0;
		for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
			double weight = 0.5/(n_aug_+lambda_aug_);
			weights_(i) = weight;
		}

		R_radar_ <<    std_radr_*std_radr_, 0, 0,
		 0, std_radphi_*std_radphi_, 0, 
		 0, 0, std_radrd_*std_radrd_;

		R_laser_ << std_laspx_*std_laspx_, 0,
		 0, std_laspy_*std_laspy_;

		previous_timestamp_ = meas_package.timestamp_;
		// done initializing, no need to predict or update
		is_initialized_ = true;
		//std::cout << "finished is_initialized =: " << is_initialized_ << endl;
		return;

	}// end initialization

	double dt =  (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;
	
	Prediction(dt);
	//std::cout << "postprediction" << endl;




	if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
		//cout << "RADAR Sensor Type" << endl;
		R_ = R_radar_;
		UpdateRadar(meas_package);
		// Radar updates
	} else {
		//cout << "LIDAR Sensor Type" << endl;
		R_ = R_laser_;
		UpdateLidar(meas_package);
	}

}// end processmeasurement


void UKF::GenerateSigmaPoints() {
	///* Sigma point spreading parameter

	//calculate square root of P
	MatrixXd A = P_.llt().matrixL();

	//set first column of sigma point matrix
	Xsig.col(0)  = x_;

	//set remaining sigma points
	for (int i = 0; i < n_x_; i++)
	{
		Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
		Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
	}
	//std::cout << "GenerateSigmaPOints x_ = " << std::endl << x_ << std::endl;
	//std::cout << "GenerateSigmaPOints lambda_ = " << std::endl << lambda_ << std::endl;
	//std::cout << "GenerateSigmaPOints n_x_ = " << std::endl << n_x_ << std::endl;
	//std::cout << "GenerateSigmaPOints A = " << std::endl << A << std::endl;
	//std::cout << "GenerateSigmaPOints Xsig = " << std::endl << Xsig << std::endl;
}
void UKF::AugmentedSigmaPoints() {




	//create augmented mean state
	x_aug.head(5) = x_;
	x_aug(5) = 0;
	x_aug(6) = 0;

	//create augmented covariance matrix
	P_aug.fill(0.0);
	P_aug.topLeftCorner(5,5) = P_;
	P_aug(5,5) = std_a_*std_a_;
	P_aug(6,6) = std_yawdd_*std_yawdd_;

	//create square root matrix
	MatrixXd L = P_aug.llt().matrixL();

	//create augmented sigma points
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i< n_aug_; i++)
	{
		Xsig_aug.col(i+1)       = x_aug + sqrt(lambda_aug_ + n_aug_) * L.col(i);
		Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug_ + n_aug_) * L.col(i);
	}

	//print result
	//std::cout << "Augmented Sigma POints = " << std::endl << Xsig_aug << std::endl;


}


void UKF::SigmaPointPrediction(double delta_t) {

	//predict sigma points
	for (int i = 0; i< 2*n_aug_+1; i++)
	{
		//extract values for better readability
		double p_x = Xsig_aug(0,i);
		double p_y = Xsig_aug(1,i);
		double v = Xsig_aug(2,i);
		double yaw = Xsig_aug(3,i);
		double yawd = Xsig_aug(4,i);
		double nu_a = Xsig_aug(5,i);
		double nu_yawdd = Xsig_aug(6,i);

		//predicted state values
		double px_p, py_p, yaw_p ;

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
			yaw_p = yaw + yawd*delta_t;
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
			yaw_p = yaw;
		}

		double v_p = v;
		double yawd_p = yawd;

		//add noise
		px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		v_p = v_p + nu_a*delta_t;
		yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred_(0,i) = px_p;
		Xsig_pred_(1,i) = py_p;
		Xsig_pred_(2,i) = v_p;
		Xsig_pred_(3,i) = yaw_p;
		Xsig_pred_(4,i) = yawd_p;
	}

	//print result
	//std::cout << "Predicted Sigma Points = " << std::endl << Xsig_pred_ << std::endl;

}

void UKF::PredictMeanAndCovariance() {

	//std::cout << "begin predictionmeancovariance " <<  std::endl;

	//predicted state mean
	x_pred_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
		x_pred_ = x_pred_ + weights_(i) * Xsig_pred_.col(i);
	}
	//std::cout << "x_pred_ predictionmeancovariance done " <<  std::endl;
	//predicted state covariance matrix
	P_pred_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

											   // state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_pred_;
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.0*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.0*M_PI;

		P_pred_ = P_pred_ + weights_(i) * x_diff * x_diff.transpose() ;
	}
	//std::cout << "p_pred_ predictionmeancovariance done " <<  std::endl;
	//print result
	//std::cout << "Predicted state" << std::endl;
	//std::cout << x_pred_ << std::endl;
	//std::cout << "Predicted covariance matrix" << std::endl;
	//std::cout << P_pred_ << std::endl;

}



/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

	//std::cout << "Begin prediction step" << endl;
	//GenerateSigmaPoints();
	AugmentedSigmaPoints();
	SigmaPointPrediction(delta_t);// Predict sigma points
	PredictMeanAndCovariance();// state and state covariance
	//std::cout << "End prediction step" << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {


	int n_z = 2;
	VectorXd z = VectorXd(n_z);
	z <<
		meas_package.raw_measurements_(0),
		meas_package.raw_measurements_(1);
	VectorXd z_pred = H_laser_ * x_pred_;
	VectorXd z_diff = z - z_pred;
	MatrixXd Ht = H_laser_.transpose();
	MatrixXd S = H_laser_ * P_pred_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_pred_ * Ht;
	MatrixXd K = PHt * Si;

	x_ = x_pred_ + (K * z_diff);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_laser_) * P_pred_;



	NIS_laser_ = z_diff.transpose()*S.inverse()*z_diff;

	//std::cout << "End Update Lidar step" << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
	//std::cout << "Begin Update Radar step" << endl;

  // Predict Radar Measurement
  //set measurement dimension, radar can measure r, phi, and r_dot
	int n_z = 3;
	VectorXd z = VectorXd(n_z);
	
	z <<
		meas_package.raw_measurements_(0),
		meas_package.raw_measurements_(1),
		meas_package.raw_measurements_(2);
	

	////set vector for weights_
	//VectorXd weights_ = VectorXd(2*n_aug_+1);
	//double weight_0 = lambda_aug_/(lambda_aug_+n_aug_);
	//weights_(0) = weight_0;
	//for (int i=1; i<2*n_aug_+1; i++) {  
	//	double weight = 0.5/(n_aug_+lambda_aug_);
	//	weights_(i) = weight;
	//}

	

	

	//transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

												// extract values for better readibility
		double p_x = Xsig_pred_(0, i);
		double p_y = Xsig_pred_(1, i);
		double v = Xsig_pred_(2, i);
		double yaw = Xsig_pred_(3, i);

		double v1 = cos(yaw)*v;
		double v2 = sin(yaw)*v;
		double row = sqrt(p_x*p_x + p_y*p_y);
		// measurement model
		Zsig_R_(0, i) = row;//r
		Zsig_R_(1, i) = 0.0;
		if (fabs(p_x) > .001 || fabs(p_y) > .001) {// avoid numerical issues if 
			Zsig_R_(1, i) = atan2(p_y, p_x);
		}//phi
		Zsig_R_(2, i) = 0.0;
		if (fabs(row) > .001) {
			Zsig_R_(2, i) = (p_x*v1 + p_y*v2) / row;   //r_dot
		}
	}

	//mean predicted measurement
	VectorXd z_pred = VectorXd(n_z);
	z_pred.fill(0.0);
	for (int i=0; i < 2*n_aug_+1; i++) {
		z_pred = z_pred + weights_(i) * Zsig_R_.col(i);
	}

	//measurement covariance matrix S
	MatrixXd S = MatrixXd(n_z,n_z);
	//std::cout << "S declared "<< std::endl;
	S.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
												//residual
		VectorXd z_diff = Zsig_R_.col(i) - z_pred;

		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		S = S + weights_(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	S = S + R_;

	//print result
	//std::cout << "z_pred radar: " << std::endl << z_pred << std::endl;
	//std::cout << "S radar: " << std::endl << S << std::endl;
	//std::cout << "S radar size: " << std::endl << S.size() << std::endl;

	//Update 


	//calculate cross correlation matrix
	Tc_R_.fill(0.0);
	for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

											   //residual
		VectorXd z_diff = Zsig_R_.col(i) - z_pred;
		//angle normalization
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

		// state difference
		VectorXd x_diff = Xsig_pred_.col(i) - x_pred_;///////////////////////////////////////////////// changed from x_
		//angle normalization
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

		Tc_R_ = Tc_R_ + weights_(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	MatrixXd K = Tc_R_ * S.inverse();
	//std::cout << "Updated K after radar: " << std::endl << K << std::endl;
	//residual
	VectorXd z_diff = z - z_pred;

	//angle normalization
	while (z_diff(1)> M_PI) z_diff(1)-=2.0*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1)+=2.0*M_PI;

	//update state mean and covariance matrix
	x_ = x_pred_ + K * z_diff;
	P_ = P_pred_ - K*S*K.transpose();


	//print result
	//std::cout << "Updated state x after radar: " << std::endl << x_ << std::endl;
	//std::cout << "Updated state covariance P after radar: " << std::endl << P_ << std::endl;
	//std::cout << "Updated state Radar" << std::endl;
	//std::cout << x_ << std::endl;
	//std::cout << "Updated covariance matrix Radar" << std::endl;
	//std::cout << P_ << std::endl;

	NIS_radar_ = z_diff.transpose()*S.inverse()*z_diff;
	// calcualate NIS
	//std::cout << "End Update Radar step" << endl;
}

