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

//Eigen::VectorXd vel=sol_darcy.head(data_d.Nx+1);

Data_Transport data_t("data.pot");
Data_Reaction  data_r("data.pot");

//Eigen::VectorXd vel=0.0*Eigen::VectorXd::Ones(data_t.Nx+1);//Per osservare il decadimento lineare
Eigen::VectorXd vel=0.0*Eigen::VectorXd::Ones(data_t.Nx+1);
Eigen::MatrixXd Ca(data_t.Nx,data_t.Nt);
Eigen::MatrixXd CaSiO3(data_t.Nx,data_t.Nt);

//std::cout<<vel.array().pow(2)<<std::endl;

Transport_system_implicit(Ca,CaSiO3,vel,data_t,data_r);

Transport_output_results_fixed_time(Ca,CaSiO3,data_t.Nx,data_t.L,data_t.Nt);

Transport_output_results_fixed_space(Ca,CaSiO3,data_t.Nx,data_t.T,data_t.Nt);

}
