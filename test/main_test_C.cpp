#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"


int main(int argc, char **argv)
{
Data data("data.pot");

auto &[domain_length, K, phi, mu, Q_in, Q_out, p_in, p_out, f, Nx, BC_in, BC_out, dt, final_time, C_in, C_out, bc_cond]=data;

std::cout<<bc_cond<<std::endl;

Eigen::VectorXd vel=2*Eigen::VectorXd::Ones(Nx+1);
double h =static_cast<double>(domain_length)/Nx;
Matrix_C C(Nx,Nx,bc_cond,phi,h,C_in);
C.set_matrix();
C.set_rhs();
C.set_BC();
C.print_m();



return 0;
} 

