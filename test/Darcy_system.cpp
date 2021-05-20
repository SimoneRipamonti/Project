#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"

int main(int argc, char **argv)
{
Data_Darcy data("data.pot");
Eigen::MatrixXd M;
Eigen::VectorXd rhs;
set_Darcy_system(data,M,rhs);
Eigen::VectorXd sol=M.fullPivLu().solve(rhs);


Darcy_output_results(sol,data.Nx,data.L);



return 0;
} 
