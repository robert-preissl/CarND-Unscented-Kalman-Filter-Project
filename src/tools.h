
#ifndef TOOLS_H_
#define TOOLS_H_


#include <vector>
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  VectorXd CalculateRMSE(const vector<VectorXd> &estimations, const vector<VectorXd> &ground_truth, int dim);

  /**
  * A helper method to convert polar coordinates into cartesian coordinates
  */
  VectorXd ConvertPolarToCartesian(const VectorXd& polar_vec);

  // constrain an angle between -pi and pi
  float constrainAngle(float x);

};

#endif /* TOOLS_H_ */
