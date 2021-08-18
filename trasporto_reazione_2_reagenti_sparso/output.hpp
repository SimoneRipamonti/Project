#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"


void output_results_fixed_time_2_reagents(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2, unsigned int Nx, double L,unsigned int Nt);

void output_results_fixed_space_2_reagents(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2,unsigned int Nx, double T,unsigned int Nt);


#endif

