#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

void set_Darcy_system(Data &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs);

void output_results(Eigen::VectorXd &value, Eigen::VectorXd &space, std::string &name);

#endif

