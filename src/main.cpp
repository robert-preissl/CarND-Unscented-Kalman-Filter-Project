
//        UKF

// state vector x = p_x, p_y, v, ψ, ψ˙
// .. v and ψ˙  are const

// transformation: x_k+1 = f(x_k + nu_k) .. non-linear

// sigma_points : representation of whole distribution. around mean and in certain relation to sigma of each state dimension
// -> throw them into f and compute mean and covariance. gives good approximation of real distribution of measurements

//

//  A) Predict

// 1. generate sigma points
//      Xk∣k=[Xk∣kX k∣k +√(λ+nx)Pk∣k  Xk∣k −√(λ+nx)Pk∣k]    // remember that Xk∣k is the first column of the Sigma matrix.

// 2. predict sigma points (15 -- see below why) - insert into process function

// 3. calculate predicted mean and covariance from predicted sigma points

// note, in process function we have noise vector. since f is non-linear, noise is applied non lin.
//  -> different to EKF. process noise is vector of 2. appended to state vector (in EKF added); plus we adjust P.
//  -> augmentation
//  -> adapt sigma point generation since now (5+2)*2 +1 points


//   B) Update

//  1. transform predicted state into meass. space;

//  2. reuse sigma points and transfer into meas space and calculate mean and covariance of predicted measurement

//  3. add meas cov noise R

//  need cross correlation between sigma points in state space and sigma points in meas space




#include <iostream>
#include "Eigen/Dense"
#include <vector>
#include "ukf.h"
#include "measurement_package.h"
#include <fstream>
#include <sstream>
#include <stdlib.h>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions1 = "Usage instructions: ";
  usage_instructions1 += argv[0];
  usage_instructions1 += " path/to/input.txt output.txt";

  string usage_instructions2 = "Usage instructions: ";
  usage_instructions2 += argv[0];
  usage_instructions2 += " path/to/input.txt output.txt skip_laser skip_radar";

  bool has_valid_args = false;

  // make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions1 << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions1 << endl;
  } else if (argc == 4) {
      cerr << "Please include an output file.\n" << usage_instructions2 << endl;
  } else if (argc == 3 || argc == 5) {
    has_valid_args = true;
  } else if (argc > 5) {
    cerr << "Too many arguments.\n" << usage_instructions1 << endl;
  }

  if (!has_valid_args) {
    exit(EXIT_FAILURE);
  }
}

void check_files(ifstream& in_file, string& in_name,
                 ofstream& out_file, string& out_name) {
  if (!in_file.is_open()) {
    cerr << "Cannot open input file: " << in_name << endl;
    exit(EXIT_FAILURE);
  }

  if (!out_file.is_open()) {
    cerr << "Cannot open output file: " << out_name << endl;
    exit(EXIT_FAILURE);
  }
}

int main(int argc, char* argv[]) {

  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  int skip_laser = 0;
  int skip_radar = 0; // by default we always use laser and rader measurements
  if(argc == 5) {
    skip_laser = atoi(argv[3]);
    skip_radar = atoi(argv[4]);
  }

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  /**********************************************
   *  Set Measurements                          *
   **********************************************/

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;
  string line;

  // prep the measurement packages (each line represents a measurement at a
  // timestamp)
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0 && skip_laser == 0) {
      // laser measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0 && skip_radar == 0) {
      // radar measurement

      // read measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float theta;
      float ro_dot;
      iss >> ro;
      iss >> theta;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, theta, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
      measurement_pack_list.push_back(meas_package);
    }

    // read ground truth data to compare later
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;

    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;

    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);

  }

  // Create a UKF instance
  UKF ukf;

  // used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  size_t number_of_measurements = measurement_pack_list.size();

  // start filtering from the second frame (the speed is unknown in the first
  // frame)
  for (size_t k = 0; k < number_of_measurements; ++k) {
    // Call the UKF-based fusion

    if(ukf.ProcessMeasurement(measurement_pack_list[k]) < 0) {
      continue; // if we got a negative responese code (e.g., 0 values), we try another measurment to initialize (only return <0 in init phase)
    }

    // output the estimation
    out_file_ << ukf.x_(0) << "\t"; // pos1 - est
    out_file_ << ukf.x_(1) << "\t"; // pos2 - est
    out_file_ << ukf.x_(2) << "\t"; // vel_abs -est
    out_file_ << ukf.x_(3) << "\t"; // yaw_angle -est
    out_file_ << ukf.x_(4) << "\t"; // yaw_rate -est

    // output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER && skip_laser == 0) {
      // output the estimation

      // p1 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";

      // p2 - meas
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR && skip_radar == 0) {
      // output the estimation in the cartesian coordinates
      float ro = measurement_pack_list[k].raw_measurements_(0);
      float phi = measurement_pack_list[k].raw_measurements_(1);
      out_file_ << ro * cos(phi) << "\t"; // p1_meas
      out_file_ << ro * sin(phi) << "\t"; // p2_meas
    }

    out_file_ << gt_pack_list[k].gt_values_(0) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(1) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(2) << "\t";
    out_file_ << gt_pack_list[k].gt_values_(3) << "\t";

    out_file_ << ukf.NIS_laser_ << "\t";
    out_file_ << ukf.NIS_radar_ << "\n";

    VectorXd positon_estimate = VectorXd(4);
    float px = ukf.x_(0);
    float py = ukf.x_(1);
    float vx = ukf.x_[2] * cos(ukf.x_[3]);
    float vy = ukf.x_[2] * sin(ukf.x_[3]);
    positon_estimate << px, py, vx, vy;

    estimations.push_back(positon_estimate);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }

  // compute the accuracy (RMSE)
  Tools tools;

  cout << " estimations.size = " << estimations.size() << " / estimations.rows() = " << estimations[0].rows() << " / estimations.cols() = " <<  estimations[0].cols() << endl;
  cout << " ground_truth.size = " << ground_truth.size() << " / ground_truth.rows() = " << ground_truth[0].rows() << " / ground_truth.cols() = " <<  ground_truth[0].cols() << endl;

  cout << "Accuracy - RMSE:" << endl << tools.CalculateRMSE(estimations, ground_truth, 4) << endl;


  // close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done!" << endl;
  return 0;
}
