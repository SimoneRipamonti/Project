#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"



//Function that saves the solution matrix (each row represent a spatial position and each column a temporal one) on a csv file 
void output_results_fixed_time(const std::string& title, const Eigen::MatrixXd &value1, unsigned int Nx, double L,unsigned int Nt);

void output_results_fixed_space(const std::string& title, const Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);


#endif

