#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using namespace std;

namespace Tools {
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations,
                         const vector<VectorXd> &ground_truth) {
    VectorXd rmse(4);
    rmse << 0, 0, 0, 0;

    // Check inputs
    if (estimations.size() != ground_truth.size()
        || estimations.size() == 0) {
      cout << "Invalid estimation or ground_truth data" << endl;
      return rmse;
    }

    // Calculate RMSE
    for (int i = 0; i < estimations.size(); ++i) {
      VectorXd res = estimations[i] - ground_truth[i];
      res = res.array() * res.array();
      rmse += res;
    }

    rmse /= estimations.size();
    rmse = rmse.array().sqrt();

    return rmse;
  }

  void polar_to_cartesian(const MeasurementPackage &measurement_pack, double &px, double &py) {
    double phi = measurement_pack.raw_measurements_[1];

    px = measurement_pack.raw_measurements_[0] * cos(phi);
    py = measurement_pack.raw_measurements_[0] * sin(phi);
  }

  VectorXd cartesian_to_polar(const Eigen::VectorXd x) {
    VectorXd z_pred(3);

    // Unpack the state vector
    double px = x(0);
    double py = x(1);
    double vx = x(2);
    double vy = x(3);

    if (fabs(px) < APPROX_ZERO) {
      px = APPROX_ZERO;
    }

    // Convert from cartesian to polar
    double px2 = px * px;
    double py2 = py * py;
    double rho = sqrt(px2 + py2);

    // Avoid division by zero
    if (fabs(rho) < APPROX_ZERO) {
      rho = APPROX_ZERO;
    }

    z_pred[0] = rho;
    z_pred[1] = atan2(py, px);
    z_pred[2] = (px * vx + py * vy) / rho;

    return z_pred;
  }

  double constrainAngle(double angle){
    angle = fmod(angle + M_PI,2*M_PI);
    if (angle < 0)
      angle += 2*M_PI;
    return angle - M_PI;
  }
}
