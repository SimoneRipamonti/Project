#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

/*!
 *Function that saves the solution vector of the Darcy system on a csv file
*/
void Darcy_output_results(Eigen::VectorXd &value,unsigned int Nx,double L);

/*!
 *Function that saves the exact pressure values at the nodes point in a csv file
*/
void pressure_exact_result(Eigen::VectorXd &value1, unsigned int Nx, double L);


/*!
 *Function that saves the exact velocity values at the nodes point in a csv file
*/
void velocity_exact_result(Eigen::VectorXd &value1, unsigned int Nx, double L);

/*!
 *Function that saves the error values in a csv file 
*/
void output_error(Eigen::VectorXd &value1, Eigen::VectorXi &Nx);


#endif

