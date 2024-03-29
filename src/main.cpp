#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include "Eigen/Dense"
#include "ukf.h"
#include "ground_truth_package.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

void check_arguments(int argc, char* argv[]) {
  string usage_instructions = "Usage instructions: ";
  usage_instructions += argv[0];
  usage_instructions += " path/to/input.txt output.txt";

  bool has_valid_args = false;

  // Make sure the user has provided input and output files
  if (argc == 1) {
    cerr << usage_instructions << endl;
  } else if (argc == 2) {
    cerr << "Please include an output file.\n" << usage_instructions << endl;
  } else if (argc == 3) {
    has_valid_args = true;
  } else if (argc > 3) {
    cerr << "Too many arguments.\n" << usage_instructions << endl;
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

  // Check files and parameters
  check_arguments(argc, argv);

  string in_file_name_ = argv[1];
  ifstream in_file_(in_file_name_.c_str(), ifstream::in);

  string out_file_name_ = argv[2];
  ofstream out_file_(out_file_name_.c_str(), ofstream::out);

  check_files(in_file_, in_file_name_, out_file_, out_file_name_);

  vector<MeasurementPackage> measurement_pack_list;
  vector<GroundTruthPackage> gt_pack_list;
  string line;

  // Prepare measurement packages for each line in the file
  while (getline(in_file_, line)) {
    string sensor_type;
    MeasurementPackage meas_package;
    GroundTruthPackage gt_package;
    istringstream iss(line);
    long long timestamp;

    // reads first element from the current line
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      // Read laser measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);

      float px, py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;

      iss >> timestamp;
      meas_package.timestamp_ = timestamp;

      measurement_pack_list.push_back(meas_package);
    } else if (sensor_type.compare("R") == 0) {
      // Read radar measurements at this timestamp
      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);

      float ro, phi, ro_dot;
      iss >> ro;
      iss >> phi;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, phi, ro_dot;

      iss >> timestamp;
      meas_package.timestamp_ = timestamp;

      measurement_pack_list.push_back(meas_package);
    }

    // Read ground truth data to compare later
    float x_gt, y_gt, vx_gt, vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;

    gt_package.gt_values_ = VectorXd(4);
    gt_package.gt_values_ << x_gt, y_gt, vx_gt, vy_gt;
    gt_pack_list.push_back(gt_package);
  }

  // Create a UKF instance
  UKF ukf(true, true);

  // Used to compute the RMSE later
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;
  size_t number_of_measurements = measurement_pack_list.size();

  // column names for output file
  out_file_ << "px" << "\t";
  out_file_ << "py" << "\t";
  out_file_ << "v" << "\t";
  out_file_ << "yaw_angle" << "\t";
  out_file_ << "yaw_rate" << "\t";
  out_file_ << "px_measured" << "\t";
  out_file_ << "py_measured" << "\t";
  out_file_ << "px_true" << "\t";
  out_file_ << "py_true" << "\t";
  out_file_ << "v_true" << "\t";
  out_file_ << "yaw_angle_true" << "\t";
  out_file_ << "yaw_rate_true" << "\t";
  out_file_ << "vx_true" << "\t";
  out_file_ << "vy_true" << "\t";
  out_file_ << "NIS" << "\n";

  // Apply the Kalman Filter to each measurement package
  for (size_t k = 0; k < number_of_measurements; ++k) {
    // Call the UKF-based fusion
    bool res = ukf.ProcessMeasurement(measurement_pack_list[k]);
    if (res==false)
      continue;

    // Output the estimation
    out_file_ << ukf.getx_()(0) << "\t"; // pos1 - est
    out_file_ << ukf.getx_()(1) << "\t"; // pos2 - est
    out_file_ << ukf.getx_()(2) << "\t"; // vel_abs - est
    out_file_ << ukf.getx_()(3) << "\t"; // yaw_angle - est
    out_file_ << ukf.getx_()(4) << "\t"; // yaw_rate - est

    // Output the measurements
    if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::LASER) {
      // output the estimation
      out_file_ << measurement_pack_list[k].raw_measurements_(0) << "\t";
      out_file_ << measurement_pack_list[k].raw_measurements_(1) << "\t";
    } else if (measurement_pack_list[k].sensor_type_ == MeasurementPackage::RADAR) {
      // output the estimation in the cartesian coordinates

      double px, py;
      Tools::polar_to_cartesian(measurement_pack_list[k], px, py);
      out_file_ << px << "\t"; // p1_meas
      out_file_ << py << "\t"; // p2_meas
    }

    // Output the ground truth packages
    double x_gt = gt_pack_list[k].gt_values_(0);
    double y_gt = gt_pack_list[k].gt_values_(1);
    double vx_gt = gt_pack_list[k].gt_values_(2);
    double vy_gt = gt_pack_list[k].gt_values_(3);
    double v_gt = sqrt(vx_gt * vx_gt + vy_gt * vy_gt);
    double yaw_gt = fabs(vx_gt) > APPROX_ZERO ? atan(vy_gt / vx_gt) : 0;
    double yaw_rate_gt = 0;

    out_file_ << x_gt << "\t";
    out_file_ << y_gt << "\t";
    out_file_ << v_gt << "\t";
    out_file_ << yaw_gt << "\t";
    out_file_ << yaw_rate_gt << "\t";
    out_file_ << vx_gt << "\t";
    out_file_ << vy_gt << "\t";

    // Output the NIS values
    out_file_ << ukf.getNIS_laser_() << "\t";
    out_file_ << ukf.getNIS_radar_() << endl;

    // Convert ukf x vector to cartesian to compare to ground truth
    VectorXd ukf_x_cartesian_ = VectorXd(4);

    float x_estimate_ = ukf.getx_()(0);
    float y_estimate_ = ukf.getx_()(1);
    float vx_estimate_ = ukf.getx_()(2) * cos(ukf.getx_()(3));
    float vy_estimate_ = ukf.getx_()(2) * sin(ukf.getx_()(3));
    
    ukf_x_cartesian_ << x_estimate_, y_estimate_, vx_estimate_, vy_estimate_;

    estimations.push_back(ukf_x_cartesian_);
    ground_truth.push_back(gt_pack_list[k].gt_values_);
  }

  // Compute the accuracy (RMSE)
  cout << endl << "Accuracy - RMSE:" << endl << Tools::CalculateRMSE(estimations, ground_truth) << endl;

  // Close files
  if (out_file_.is_open()) {
    out_file_.close();
  }

  if (in_file_.is_open()) {
    in_file_.close();
  }

  cout << "Done processing " << in_file_name_ << endl;
  return 0;
}
