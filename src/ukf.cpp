#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const double INITIAL_VAR = 0.5;      // large variance
#define COL_PHI 1
#define COL_PSI 3
#define MAX_DT_DIFF 0.1
#define DT_STEP MAX_DT_DIFF/2

// https://discussions.udacity.com/t/numerical-instability-of-the-implementation/230449/11

UKF::UKF(bool enable_radar, bool enable_lidar) {
  use_laser_ = enable_radar;
  use_radar_ = enable_lidar;

  is_initialized_ = false;

  n_x_ = 5;
  n_aug_ = 7;
  n_sigma = 2 * n_aug_ + 1;
  n_z_radar = 3;
  n_z_lidar = 2;

  x_ = VectorXd(n_x_);
  P_ = MatrixXd(n_x_, n_x_);
  P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

  lambda_ = 3 - n_aug_;

  // Process noise
  std_a_ = 2;         // good: 2
  std_yawdd_ = 0.5;   // good: 0.5

  // Measurement noise - lidar
  std_laspx_ = 0.15;
  std_laspy_ = 0.15;

  // Measurement noise - radar
  std_radr_ = 0.3;
  std_radphi_ = 0.03;     // 0.0175
  std_radrd_ = 0.3;       // 0.1

  H_laser_ = MatrixXd(n_z_lidar, n_x_);
  H_laser_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  R_laser_ = MatrixXd(n_z_lidar, n_z_lidar);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  R_radar_ = MatrixXd(n_z_radar,n_z_radar);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;

  Q_ = MatrixXd::Zero(n_z_lidar, n_z_lidar);
  Q_ << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;

  // Init weight vector
  weights_ = VectorXd(n_sigma);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<n_sigma; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

}

UKF::~UKF() {}

// Getters
Eigen::VectorXd UKF::getx_() {
  return x_;
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  // Perform one-off initialization using the first measurement
  if (!is_initialized_) {
    double px;
    double py;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar measurement from polar to cartesian coordinates
      Tools::polar_to_cartesian(meas_package, px, py);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // No conversion required for lidar measurements
      px = meas_package.raw_measurements_[0];
      py = meas_package.raw_measurements_[1];
    }

    // Set initial state x to the first measurement's position and zero velocity
    x_ << px, py, 0, 0, 0;

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  // 1. Prediction step
  // Compute the time elapsed in seconds between the current and previous measurements
  float dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // Predict state vector x and state covariance matrix P, watch out for large dt values
  while (dt > MAX_DT_DIFF)
  {
    Prediction(DT_STEP);
    dt -= DT_STEP;
  }
  MatrixXd Xsig_pred = Prediction(dt);

  // 2. Measurement update/correction step
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package, Xsig_pred);
  } else {
    UpdateLidar(meas_package);
  }
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
MatrixXd UKF::Prediction(double delta_t) {
  MatrixXd Xsig_aug, Xsig_pred;

  //  Xsig = GenerateSigmaPoints();
  Xsig_aug = AugmentedSigmaPoints();
  Xsig_pred = SigmaPointPrediction(Xsig_aug, delta_t);
  PredictMeanAndCovariance(Xsig_pred);

  return Xsig_pred;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  You'll also need to calculate the lidar NIS.
  */
  VectorXd z_pred;
  MatrixXd S, Zsig;

//  PredictLaserMeasurement(Xsig_pred, &z_pred, &S, &Zsig);
  UpdateStateLaser(meas_package.raw_measurements_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package, MatrixXd Xsig_pred) {
  /**
  TODO You'll also need to calculate the radar NIS.
  */
  VectorXd z_pred;
  MatrixXd S, Zsig;

  PredictRadarMeasurement(Xsig_pred, &z_pred, &S, &Zsig);
  UpdateStateRadar(Xsig_pred, Zsig, z_pred, S, meas_package.raw_measurements_);
}

MatrixXd UKF::AugmentedSigmaPoints() {
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sigma);

  // Create augmented mean state
  x_aug << x_, 0, 0;

  P_aug.fill(0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug.bottomRightCorner(2,2) = Q_;

  // Create square root matrix
  MatrixXd A = P_aug.llt().matrixL();

  // Create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 1; i <= n_aug_; i++)
  {
    Xsig_aug.col(i)        = x_aug + sqrt(lambda_+n_aug_) * A.col(i-1);
    Xsig_aug.col(i+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * A.col(i-1);
  }

  return Xsig_aug;
}

MatrixXd UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sigma);
  Xsig_pred.fill(0);

  //predict sigma points
  for (int i=0;i < n_sigma;i++) {
    double vk = Xsig_aug(2,i);
    double psi = Xsig_aug(3,i);
    double psidot = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_psidotdot = Xsig_aug(6,i);

    if (fabs(psidot) > 0.001) {
      Xsig_pred(0, i) = Xsig_aug(0, i) + (vk / psidot) * (sin(psi + psidot * delta_t) - sin(psi)) +
                        0.5 * delta_t * delta_t * cos(psi) * nu_a;
      Xsig_pred(1, i) = Xsig_aug(1, i) + (vk / psidot) * (-cos(psi + psidot * delta_t) + cos(psi)) +
                        0.5 * delta_t * delta_t * sin(psi) * nu_a;
    }
    else {
      Xsig_pred(0, i) = Xsig_aug(0, i) + (vk * cos(psi)*delta_t) + 0.5 * delta_t * delta_t * cos(psi) * nu_a;
      Xsig_pred(1, i) = Xsig_aug(1, i) + (vk * sin(psi)*delta_t) + 0.5 * delta_t * delta_t * sin(psi) * nu_a;
    }
    Xsig_pred(2,i) = Xsig_aug(2, i) + 0 + delta_t*nu_a;
    Xsig_pred(3,i) = Xsig_aug(3, i) + psidot*delta_t + 0.5*delta_t*delta_t*nu_psidotdot;
    Xsig_pred(4,i) = Xsig_aug(4, i) + 0 + delta_t*nu_psidotdot;
  }

  return Xsig_pred;
}

void UKF::PredictMeanAndCovariance(MatrixXd Xsig_pred) {
  // Predict state mean
  x_.fill(0);
  for(int i=0;i<n_sigma;i++) {
    x_+= weights_(i) * Xsig_pred.col(i);
  }

  // Predict state covariance matrix
  P_.fill(0);
  MatrixXd prod;
  for(int i=0;i<n_sigma;i++) {
    prod = Xsig_pred.col(i) - x_;

    // Normalise psi to be between -pi and pi
    while (prod(COL_PSI)> M_PI) prod(COL_PSI)-=2.*M_PI;
    while (prod(COL_PSI)< -M_PI) prod(COL_PSI)+=2.*M_PI;

    P_ += weights_(i) * prod * prod.transpose();
  }
}

void UKF::PredictRadarMeasurement(MatrixXd Xsig_pred, VectorXd* z_pred_out, MatrixXd* S_out, MatrixXd* Zsig_out) {
  MatrixXd Zsig = MatrixXd(n_z_radar, n_sigma);     // Sigma points in measurement space
  VectorXd z_pred = VectorXd(n_z_radar);            // Mean predicted measurement
  MatrixXd S = MatrixXd(n_z_radar, n_z_radar);      // Measurement covariance matrix S

  // Transform sigma points into measurement space
  for (int i;i<n_sigma;i++) {
    double px = Xsig_pred(0,i);
    double py = Xsig_pred(1,i);
    double v = Xsig_pred(2,i);
    double yaw = Xsig_pred(3, i);

    double vx = cos(yaw)*v;
    double vy = sin(yaw)*v;
    double rho = sqrt(px*px + py*py);

    Zsig(0, i) = rho;
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * vx + py * vy) / rho;
  }

  // Calculate mean predicted measurement
  z_pred.fill(0);
  for(int i=0;i<n_sigma;i++) {
    z_pred+= weights_(i) * Zsig.col(i);
  }

  // Calculate measurement covariance matrix S
  S.fill(0);
  MatrixXd prod;
  for(int i=0;i<n_sigma;i++) {
    prod = Zsig.col(i) - z_pred;
    // normalise phi to be between -pi and pi
    while (prod(COL_PHI)> M_PI) prod(COL_PHI)-=2.*M_PI;
    while (prod(COL_PHI)<-M_PI) prod(COL_PHI)+=2.*M_PI;

    S += weights_(i) * prod * prod.transpose();
  }
  S = S + R_radar_;

  // Write result
  *z_pred_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

void UKF::UpdateStateRadar(MatrixXd Xsig_pred, MatrixXd Zsig, VectorXd z_pred, MatrixXd S, VectorXd z) {
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar);

  // Calculate cross correlation matrix
  Tc.fill(0);
  for (int i=0;i<n_sigma;i++) {
    MatrixXd first = Xsig_pred.col(i) - x_;
    while (first(COL_PSI)> M_PI) first(COL_PSI)-=2.*M_PI;
    while (first(COL_PSI)<-M_PI) first(COL_PSI)+=2.*M_PI;

    MatrixXd second = Zsig.col(i) - z_pred;
    while (second(COL_PHI)> M_PI) second(COL_PHI)-=2.*M_PI;
    while (second(COL_PHI)<-M_PI) second(COL_PHI)+=2.*M_PI;

    Tc += weights_(i) * first * second.transpose();
  }
  // Calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  // Update state mean and covariance matrix
  VectorXd diff = z - z_pred;

  // Normalize phi
  while (diff(COL_PHI)> M_PI) diff(COL_PHI)-=2.*M_PI;
  while (diff(COL_PHI)<-M_PI) diff(COL_PHI)+=2.*M_PI;

  x_ = x_ + K * diff;
  P_ = P_ - K * S * K.transpose();
}

void UKF::UpdateStateLaser(VectorXd z) {
  // Perform various matrix calculations
  VectorXd z_pred = H_laser_ * x_;

  VectorXd y = z - z_pred;
  MatrixXd Ht = H_laser_.transpose();
  MatrixXd S = H_laser_ * P_ * Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  // New estimates for x and P
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_laser_) * P_;
}
