#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

//Function that saves the solution vector of the Darcy system on a csv file
void Darcy_output_results(Eigen::VectorXd &value,unsigned int Nx,double L);

//Function that saves the solution matrix (each row represent a spatial position and each column a temporal one) on a csv file 
void Transport_output_results_fixed_time(Eigen::MatrixXd &value1, unsigned int Nx, double L,unsigned int Nt);

void Transport_output_results_fixed_space(Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);

void Transport_output_results_fixed_space2(Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);

void Transport_output_results_fixed_space3(Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);

void Transport_output_results_fixed_space4(Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);

void Transport_output_results_fixed_space5(Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);


#endif

