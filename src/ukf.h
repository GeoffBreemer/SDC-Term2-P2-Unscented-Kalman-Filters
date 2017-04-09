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
private:
  bool is_initialized_;
  long long previous_timestamp_;

  bool use_laser_;        // if this is false, laser measurements will be ignored (except for init)
  bool use_radar_;        // if this is false, radar measurements will be ignored (except for init)

  int n_x_;               // State dimension
  int n_aug_;             // Augmented state dimension
  int n_sigma;
  int n_z_radar;          // Radar measurement dimension: r, phi, and r_dot
  int n_z_laser;          // Laser measurement dimension: px and py

  VectorXd x_;            // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  MatrixXd P_;            // state covariance matrix
  MatrixXd Q_;            // augmented covariance matrix
  double lambda_;         // Sigma point spreading parameter
  double NIS_radar_;      // the current NIS for radar
  double NIS_laser_;      // the current NIS for laser

  double std_a_;          // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_yawdd_;      // Process noise standard deviation yaw acceleration in rad/s^2

  double std_laspx_;      // Laser measurement noise standard deviation position1 in m
  double std_laspy_;      // Laser measurement noise standard deviation position2 in m

  double std_radr_;       // Radar measurement noise standard deviation radius in m
  double std_radphi_;     // Radar measurement noise standard deviation angle in rad
  double std_radrd_ ;     // Radar measurement noise standard deviation radius change in m/s

  MatrixXd R_radar_;      // Measurement covariance matrix - radar measurement noise
  MatrixXd R_laser_;      // Measurement covariance matrix - laser measurement noise
  MatrixXd H_laser_;      // Measurement function H matrix

  VectorXd w_;            // Sigma points weights

public:

  UKF(bool enable_radar = true, bool enable_lidar = true);
  virtual ~UKF();

  // Getters
  VectorXd getx_() {return x_;};
  double getNIS_radar_() {return NIS_radar_;};
  double getNIS_laser_() {return NIS_laser_;};

  // Run the Predict -> Measurement Update process for one measurement
  bool ProcessMeasurement(MeasurementPackage meas_package);

  // Predicts sigma points, the state, and the state covariance
  MatrixXd Prediction(double delta_t);

  // Generates augmented sigma points
  MatrixXd AugmentedSigmaPoints();

  // Predicts sigma points using the process model
  MatrixXd SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t);

  // Predicts state mean and covariance
  void PredictMeanAndCovariance(MatrixXd Xsig_pred);

  // Predicts radar measurement mean and covariance using the measurement model
  void PredictRadarMeasurement(MatrixXd Xsig_pred, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out);

  // Updates the state and the state covariance matrix using a radar measurement, computes NIS
  void UpdateRadar(MatrixXd Xsig_pred, MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z);

  // Updates the state and the state covariance matrix using a laser measurement, computes NIS
  void UpdateLaser(VectorXd z);
};

#endif /* UKF_H */
