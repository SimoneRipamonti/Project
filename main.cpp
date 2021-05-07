#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "functions.hpp"


int main(int argc, char **argv)
{
Data data("data.pot");
Eigen::MatrixXd M;
Eigen::VectorXd rhs;
set_Darcy_system(data,M,rhs);
Eigen::VectorXd sol=M.fullPivLu().solve(rhs);


output_results(sol,data.Nx,data.domain_length);



return 0;
} 
