#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

//Function that saves the solution vector of the Darcy system on a csv file
void Darcy_output_results(Eigen::VectorXd &value,unsigned int Nx,double L);

//Function that saves the solution matrix (each row represent a spatial position and each column a temporal one) on a csv file 
void Transport_output_results(Eigen::MatrixXd &value,unsigned int Nx, double L,unsigned int Nt);

#endif

