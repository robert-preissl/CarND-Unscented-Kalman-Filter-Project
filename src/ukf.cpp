#include <iostream>
#include "ukf.h"
#include "tools.h"

using namespace std;

/**
 * Initializes the Unscented Kalman Filter - UKF
 */
UKF::UKF() {

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  //Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.05;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.1;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.01;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.1;

  is_initialized_ = false;

  // set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //define spreading parameter
  lambda_ = 3 - n_aug_;

  n_z_radar_ = 3;
  n_z_lidar_ = 2;

  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2*n_aug_+1);

  tools = new Tools();
}

int UKF::Initialize(const MeasurementPackage& meas_package) {

  /**
    * Initialize the state ukf.x_ with the first measurement.
    * Create the covariance matrix.
    * !convert radar from polar to cartesian coordinates.
  */

  // first measurement
  cout << "UKF initialization " << endl;
  VectorXd x_init = VectorXd(n_x_);

  // state covariance matrix P
  MatrixXd P_init = MatrixXd(n_x_, n_x_);
  // initial covariance matrix
  P_init << 1, 0, 0, 0, 0,
            0, 1, 0, 0, 0,
            0, 0, 1, 0, 0,
            0, 0, 0, 1, 0,
            0, 0, 0, 0, 1;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    /**
    Convert radar from polar to cartesian coordinates to set init state.
    */
    VectorXd position_polar_coords(2);
    VectorXd position_cartesian_coords(2);
    // 0 .. radial distance from origin (ro) || 1 .. angle (theta)
    position_polar_coords << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];

    position_cartesian_coords = tools->ConvertPolarToCartesian(position_polar_coords);

    // state vector x = p_x, p_y, v, ψ, ψ˙
    x_init << position_cartesian_coords(0), position_cartesian_coords(1), 0.0, 0.0, 0.0;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    x_init << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0.0, 0.0, 0.0;
  }

  if (x_init(0) == 0 || x_init(1) == 0) {
    return -1;
  }

  // noise covariance matrix for Radar
  MatrixXd R_radar_init = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_init << std_radr_ * std_radr_, 0, 0,
                  0, std_radphi_ * std_radphi_, 0,
                  0, 0, std_radrd_ * std_radrd_;

  // noise covariance matrix for Laser
  MatrixXd R_laser_init = MatrixXd(n_z_lidar_, n_z_lidar_);
  R_laser_init << std_laspx_*std_laspx_,0,
                  0,std_laspy_*std_laspy_;

  // set initial weights
  VectorXd weights_init = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_ + n_aug_);
  weights_init(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; ++i) {  //2n+1 weights
    weights_init(i) = 0.5/(n_aug_+lambda_);
  }

  // call the Kalman filter init function to set initial vectors and matrices
  AssignInitValues(x_init, P_init, R_radar_init, R_laser_init, weights_init);

  previous_timestamp_ = meas_package.timestamp_;

  // done initializing, no need to predict or update
  is_initialized_ = true;

  cout << "UKF initialization done " << endl;

  return 0;
}

// assign initial values
void UKF::AssignInitValues(VectorXd x_in, MatrixXd P_in, MatrixXd R_radar_init, MatrixXd R_laser_init, VectorXd weights_init) {
  x_       = x_in;
  P_       = P_in;
  R_radar_ = R_radar_init;
  R_laser_ = R_laser_init;
  weights_ = weights_init;
}

UKF::~UKF() {
  delete tools;
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
int UKF::ProcessMeasurement(const MeasurementPackage& meas_package) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
   if (!is_initialized_) {
     return Initialize(meas_package);
   }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
   delta_t_ = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;
   previous_timestamp_ = meas_package.timestamp_;
   printf("delta_t_ = %.17g \n", delta_t_);

   // for large delta_t_ values call Predict repeatedly
   //  (tip taken from Slack SDCN channel)
   while (delta_t_ > 0.1) {
       double tmp = delta_t_;
       delta_t_ = 0.05;
       Predict();
       delta_t_ = tmp;
       delta_t_ -= 0.05;
   }

   Predict();


  /*****************************************************************************
   *  Update
   ****************************************************************************/
   Update(meas_package);


   // print the output
   std::cout << endl;
   std::cout << endl;

   cout << "x_ = " << x_ << endl;
   cout << "P_ = " << P_ << endl;

   return 0;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 */
void UKF::Predict() {

  // 1. create (augmented) sigma points
  MatrixXd Xsig_aug = AugmentedSigmaPoints();

  // 2. sigma point Prediction
  SigmaPointPrediction(Xsig_aug);

  // 3. calculate predicted state & covariance
  PredictMeanAndCovariance();

}

void UKF::Update(const MeasurementPackage& meas_package) {
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    UpdateLidar(meas_package);
  }
}

// compute cross correlation matrix
MatrixXd UKF::CrossCorrelationMatrix(int n_z, const MatrixXd& Zsig, const VectorXd& z_pred, bool angle_normalize) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    z_diff(1) = tools->constrainAngle(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    //angle normalization
    x_diff(3) = tools->constrainAngle(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  return Tc;
}

// compute the measurement covariance matrix
MatrixXd UKF::MeasurementCovarianceMatrix(int n_z, const MatrixXd& Zsig, const VectorXd& z_pred, bool angle_normalize) {
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
      //residual
      VectorXd z_diff = Zsig.col(i) - z_pred;
      //angle normalization
      z_diff(1) = tools->constrainAngle(z_diff(1));

      S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  return S;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(const MeasurementPackage& meas_package) {

  MatrixXd Zsig = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);
  VectorXd z_pred = VectorXd(n_z_lidar_);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }

  z_pred = Zsig * weights_;

  //measurement covariance matrix S
  MatrixXd S = MeasurementCovarianceMatrix(n_z_lidar_, Zsig, z_pred, false);

  S = S + R_laser_;

  //get matrix for cross correlation Tc
  MatrixXd Tc = CrossCorrelationMatrix(n_z_lidar_, Zsig, z_pred, false);

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = VectorXd(n_z_lidar_);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];
  VectorXd z_diff = z - z_pred;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;

  //print result
  std::cout << endl;
  std::cout << "Lidar -- Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Lidar -- Updated state covariance P: " << std::endl << P_ << std::endl;

}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(const MeasurementPackage& meas_package) {

  //create matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      // extract values for better readibility
      double p_x = Xsig_pred_(0, i);
      double p_y = Xsig_pred_(1, i);
      double v   = Xsig_pred_(2, i);
      double yaw = Xsig_pred_(3, i);

      double v_y = cos(yaw) * v;
      double v_x = sin(yaw) * v;

      Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);
      Zsig(1, i) = atan2(p_y, p_x);
      Zsig(2, i) = (p_x * v_y + p_y * v_x) / Zsig(0, i);
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);

  z_pred.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MeasurementCovarianceMatrix(n_z_radar_, Zsig, z_pred, true);

  //add measurement noise covariance matrix
  S = S + R_radar_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = CrossCorrelationMatrix(n_z_radar_, Zsig, z_pred, true);

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z = VectorXd(n_z_radar_);
  z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];
  VectorXd z_diff = z - z_pred;

  //angle normalization
  z_diff(1) = tools->constrainAngle(z_diff(1));

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();

  NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;

  //print result
  std::cout << endl;
  std::cout << "Radar -- Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Radar -- Updated state covariance P: " << std::endl << P_ << std::endl;
}


// compute augmented sigma points
MatrixXd UKF::AugmentedSigmaPoints() {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

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
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  //print result
  std::cout << endl;
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //return result
  return Xsig_aug;

}

// predict sigma points (insert augmented sigma points into state function)
void UKF::SigmaPointPrediction(const MatrixXd& Xsig_aug) {

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
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t_) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t_) );
    }
    else {
        px_p = p_x + v*delta_t_*cos(yaw);
        py_p = p_y + v*delta_t_*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t_;
    double yawd_p = yawd;

    //add noise
    double delta_t_sq = delta_t_ * delta_t_;
    px_p = px_p + 0.5*nu_a*delta_t_sq * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t_sq * sin(yaw);
    v_p = v_p + nu_a*delta_t_;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t_sq;
    yawd_p = yawd_p + nu_yawdd*delta_t_;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //print result
  std::cout << endl;
  std::cout << "Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl;
}

// calculate the predicted state covariance matrix and state x
void UKF::PredictMeanAndCovariance() {

  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  //iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  //iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = tools->constrainAngle(x_diff(3));

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }

  //print result
  std::cout << endl;
  std::cout << "Predicted state" << std::endl;
  std::cout << x_ << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P_ << std::endl;
}
