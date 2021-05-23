#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"
#include <fstream>


int main(int argc, char **argv)
{
Data_Darcy data_d("data.pot");

Eigen::MatrixXd M;
Eigen::VectorXd rhs;
set_Darcy_system(data_d,M,rhs);
Eigen::VectorXd sol_darcy=M.fullPivLu().solve(rhs);

Darcy_output_results(sol_darcy,data_d.Nx,data_d.L);

Eigen::VectorXd vel=sol_darcy.head(data_d.Nx+1);

Data_Transport data_t("data.pot");

//Eigen::VectorXd vel=1*Eigen::VectorXd::Ones(data_t.Nx+1);
Eigen::MatrixXd solution(data_t.Nx,data_t.Nt);

Transport_system_implicit(solution,vel,data_t);

Transport_output_results(solution,data_t.Nx,data_t.L,data_t.Nt);
}
