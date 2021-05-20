#ifndef SYSTEMS_HH
#define SYSTEMS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

void set_Darcy_system(Data_Darcy &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs);

void Transport_system_esplicit(Eigen::MatrixXd &solution,Eigen::VectorXd &velocity,Data_Transport &data);

#endif

