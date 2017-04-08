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
  bool use_laser_;        // if this is false, laser measurements will be ignored (except for init)
  bool use_radar_;        // if this is false, radar measurements will be ignored (except for init)

  int n_x_;               // State dimension
  int n_aug_;             // Augmented state dimension

  VectorXd x_;            // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  MatrixXd P_;            // state covariance matrix
  double lambda_;         // Sigma point spreading parameter

  double std_a_;          // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_yawdd_;      // Process noise standard deviation yaw acceleration in rad/s^2

  double std_laspx_;      // Laser measurement noise standard deviation position1 in m
  double std_laspy_;      // Laser measurement noise standard deviation position2 in m

  double std_radr_;       // Radar measurement noise standard deviation radius in m
  double std_radphi_;     // Radar measurement noise standard deviation angle in rad
  double std_radrd_ ;     // Radar measurement noise standard deviation radius change in m/s

  MatrixXd Xsig_pred_;    // predicted sigma points matrix

public:

  long long time_us_;     // time when the state is true, in us

  VectorXd weights_;      // Weights of sigma points

  double NIS_radar_;      // the current NIS for radar
  double NIS_laser_;      // the current NIS for laser

  UKF(bool enable_radar = true, bool enable_lidar = true);
  virtual ~UKF();

  // Getters
  Eigen::VectorXd getx_();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

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

  //Step 1
  void GenerateSigmaPoints(MatrixXd* Xsig_out);
  void AugmentedSigmaPoints(MatrixXd* Xsig_out);

  // Step 2
  void SigmaPointPrediction(MatrixXd* Xsig_out);

  // Step 3
  void PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred);

  // Step 4
  void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out);
  void PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out);

  // Step 5
  void UpdateState(VectorXd* x_out, MatrixXd* P_out);
};

#endif /* UKF_H */
