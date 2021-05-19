#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "functions.hpp"


int main(int argc, char **argv)
{
Data data("data.pot");

auto &[domain_length, K, phi, mu, Q_in, Q_out, p_in, p_out, f, Nx, BC_in, BC_out, dt, final_time, C_in, C_out, bc_cond]=data;

std::cout<<bc_cond<<std::endl;

Eigen::VectorXd vel=2*Eigen::VectorXd::Ones(Nx+1);
double h =static_cast<double>(domain_length)/Nx;
Matrix_F_piu F_p(Nx,Nx,bc_cond,vel);
Matrix_F_meno F_m(Nx,Nx,bc_cond,vel);
F_p.set_matrix();
F_p.set_rhs();
F_p.set_BC();
F_p.print_m();
std::cout<<" "<<std::endl;
F_m.set_matrix();
F_m.set_rhs();
F_m.set_BC();
F_m.print_m();

return 0;
} 
