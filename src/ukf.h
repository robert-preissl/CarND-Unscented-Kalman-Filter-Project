
#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include "ground_truth_package.h"
#include <vector>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  // matrix with predicted sigma points as columns
  MatrixXd Xsig_pred_;

  // noise covariance matrix for Radar & Laser
  MatrixXd R_radar_;
  MatrixXd R_laser_;

  ///* time when the state is true, in us
  long time_us_;

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

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  int n_z_radar_;
  int n_z_lidar_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  // previous timestamp
  long previous_timestamp_;

  double delta_t_;

  Tools* tools;

  /**
   * Constructor
   */
  UKF();

  /**
   * Initialization
   */
  int Initialize(const MeasurementPackage& meas_package);

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * Assign initial values
   */
  void AssignInitValues(VectorXd x_in, MatrixXd P_in, MatrixXd R_radar_init, MatrixXd R_laser_init, VectorXd weights_init);

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   * @param gt_package The ground truth of the state x at measurement time
   */
  int ProcessMeasurement(const MeasurementPackage& meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Predict();

  void Update(const MeasurementPackage& meas_package);

  // compute augmented sigma points
  MatrixXd AugmentedSigmaPoints();

  // predict sigma points
  void SigmaPointPrediction(const MatrixXd& Xsig_aug);

  // calculate the predicted state covariance matrix and state x
  void PredictMeanAndCovariance();

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(const MeasurementPackage& meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(const MeasurementPackage& meas_package);

  // compute cross correlation matrix
  MatrixXd CrossCorrelationMatrix(int n_z, const MatrixXd& Zsig, const VectorXd& z_pred, bool angle_normalize);

  // compute the measurement covariance matrix
  MatrixXd MeasurementCovarianceMatrix(int n_z, const MatrixXd& Zsig, const VectorXd& z_pred, bool angle_normalize);
};

#endif /* UKF_H */
