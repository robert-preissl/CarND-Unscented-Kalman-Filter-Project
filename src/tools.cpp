#include <iostream>
#include "tools.h"
#include <cmath>

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		                   const vector<VectorXd> &ground_truth, int dim){

	VectorXd rmse(dim);
	if (dim == 2)
		rmse << 0,0;
	else if (dim == 4)
		rmse << 0,0,0,0;
	else
		return rmse;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

VectorXd Tools::ConvertPolarToCartesian(const VectorXd& polar_vec) {
	VectorXd cartesian_vec(2);
	// http://keisan.casio.com/exec/system/1223527679
	cartesian_vec << polar_vec(0) * cos(polar_vec(1)), polar_vec(0) * sin(polar_vec(1));
	return cartesian_vec;
}

// constrain an angle between -pi and pi
// http://stackoverflow.com/questions/11498169/dealing-with-angle-wrap-in-c-code
float Tools::constrainAngle(float x) {
    x = fmod(x + M_PI, 2* M_PI);
    if (x < 0)
        x += 2*M_PI;

    return x - M_PI;
}
