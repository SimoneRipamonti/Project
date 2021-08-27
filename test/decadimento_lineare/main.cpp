//#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"
#include <fstream>

#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


int main(int argc, char **argv)
{
/*Data_Darcy data_d("data.pot");

Eigen::MatrixXd M;
Eigen::VectorXd rhs;
set_Darcy_system(data_d,M,rhs);
Eigen::VectorXd sol_darcy=M.fullPivLu().solve(rhs);

Darcy_output_results(sol_darcy,data_d.Nx,data_d.L);
*/
//Eigen::VectorXd vel=sol_darcy.head(data_d.Nx+1);

Data_Transport data_transport("data.pot");
Data_linear_decay data_linear_decay("data.pot");

//Eigen::VectorXd vel=0.0*Eigen::VectorXd::Ones(data_t.Nx+1);//Per osservare il decadimento lineare
Eigen::VectorXd vel=1.0*Eigen::VectorXd::Ones(data_transport.Nx+1);
Eigen::MatrixXd Ca(data_transport.Nx,data_transport.Nt);
//Eigen::MatrixXd CaSiO3(data_transport.Nx,data_transport.Nt);

//std::cout<<vel.array().pow(2)<<std::endl;

Transport_system_implicit(Ca,vel,data_transport,data_linear_decay);

output_results_fixed_space("Ca_fixed_space.csv", Ca, data_transport.Nx, data_transport.T, data_transport.Nt);

output_results_fixed_time("Ca_fixed_time.csv", Ca, data_transport.Nx, data_transport.L, data_transport.Nt);
}
