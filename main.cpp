#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"

void set_Darcy_system(Data &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs)
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

Eigen::MatrixXd A_=A.get_matrix();
Eigen::MatrixXd B_=B.get_matrix();
Eigen::MatrixXd B_T=(B.get_matrix()).transpose();

Eigen::MatrixXd C(A_.rows(),A_.cols()+B_T.cols());
C<<A_,B_T;

Eigen::MatrixXd D= Eigen::MatrixXd::Zero(B_.rows(),B_T.cols());
Eigen::MatrixXd E(B_.rows(),B_.cols()+D.cols());
E<<-B_,D;

Eigen::MatrixXd M(C.rows()+E.rows(),C.cols());
M<<C,E;

Eigen::VectorXd v1=A.get_rhs();
Eigen::VectorXd v2=B.get_rhs();
Eigen::VectorXd v3(v1.size()+v2.size());
v3<<v1,v2;

rhs=v3;
Matrix=M;

}

int main(int argc, char **argv)
{
Data data("data.pot");
Eigen::MatrixXd M;
Eigen::VectorXd rhs;
set_Darcy_system(data,M,rhs);

Eigen::MatrixXd sol=M.fullPivLu().solve(rhs);
std::cout<<sol<<std::endl;




/*Eigen::MatrixXd A=Eigen::MatrixXd::Zero(3,3);
Eigen::MatrixXd B=Eigen::MatrixXd::Zero(3,2);
Eigen::MatrixXd C(A.rows(),A.cols()+B.cols());
C<<A,B;
Eigen::MatrixXd D=Eigen::MatrixXd::Zero(1,5);
Eigen::MatrixXd E(C.rows()+D.rows(),D.cols());
E<<C,D;
std::cout<<E<<std::endl;

Eigen::VectorXd v(3);
v<<1,2,3;
Eigen::VectorXd v2(5);
v2<<4,5,6,7,8;
Eigen::VectorXd v3(v.size()+v2.size());
v3<<v,v2;
std::cout<<v3<<std::endl;
std::cout<<"Hello world!"<<std::endl;
std::cout<<"This is a new added line"<<std::endl;*/

return 0;
} 
