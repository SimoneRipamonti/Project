#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

/*!
 *Function that saves the solution in a .csv file the following way: each row is a spatial position, each column is a time instant 
*/
void output_results_fixed_time(const std::string& title, const Eigen::MatrixXd &value1, double L, unsigned int Nx, double T,unsigned int Nt);


/*!
 *Function that saves the solution in a .csv file the following way: each row is a time instant, each column is a spatial position 
*/
void output_results_fixed_space(const std::string& title, const Eigen::MatrixXd &value1, double L, unsigned int Nx, double T,unsigned int Nt);





#endif

