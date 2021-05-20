#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

void Darcy_output_results(Eigen::VectorXd &value,unsigned int Nx,double L);

void Transport_output_results(Eigen::MatrixXd &value,unsigned int Nx, double L,unsigned int Nt);

#endif

