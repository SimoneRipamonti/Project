#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"

Eigen::MatrixXd set_Darcy_system(Data &data)
{
auto &[domain_length, K, phi, mu, Q_in, Q_out, p_in, p_out, f, Nx, BC_in, BC_out, dt, final_time, C_in, C_out, bc_cond]=data;

double h =static_cast<double>(domain_length)/Nx;

Matrix_A A(Nx,Nx,K,h,mu,BC_in,BC_out,Q_in,p_out);
A.set_matrix();
A.set_BC();
A.set_rhs();

Matrix_B B(Nx-1,Nx,BC_in,BC_out,f,h);
B.set_matrix();
B.set_BC();
B.set_rhs();


return A.get_matrix();

}

int main(int argc, char **argv)
{
Eigen::VectorXd v(3);
v << 1, 2, 3;
std::cout<<v.tail(1)<<std::endl;
std::cout<<"Hello world!"<<std::endl;
std::cout<<"This is a new added line"<<std::endl;

return 0;
} 
