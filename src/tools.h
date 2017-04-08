#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"
#include "measurement_package.h"

const double APPROX_ZERO = 0.0001;   // avoid numerical stability issues: do not compare to zero

namespace Tools {
  //  Calculate RMSE
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  // Convert measurement from polar to cartesian coordinates
  void polar_to_cartesian(const MeasurementPackage& measurement_pack, double& px, double& py);

  // Convert state from cartesian to polar coordinates
  Eigen::VectorXd cartesian_to_polar(const Eigen::VectorXd x);
};

#endif /* TOOLS_H_ */
