#include "systems.hpp"
#include <iostream>
#include <vector>

void set_Darcy_system(Data_Darcy &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs)
{
//All the data that are needed to define the Darcy System are extracted from the data structure
auto &[L, K, phi, mu, Q_in, Q_out, p_in, p_out, f, Nx, BC_in, BC_out]=data;

//Computation of the spatial step
double h =static_cast<double>(L)/Nx;

//Assembling of the Mass Matrix of the first equation with its Boundary Condition
Matrix_A A(Nx+1,Nx+1);
A.assemble_matrix(K,h,mu,BC_in,BC_out,p_in,Q_out);

//Assembling of the B Matrix for the first equation of the system and assembling of the -B^T for the second equation of the system
//For B we have to separate all the steps (we can't use the assemble_matrix() function since for -B^T the modification due to the BC condition have not to be set)
Matrix_B B(Nx+1,Nx);//è diversa dalle note, la definisco come la B^T delle note
B.set_data(BC_in,BC_out,f,h);
B.define_matrix();
Eigen::MatrixXd B_T=(B.get_matrix()).transpose();
B.set_BC();
B.set_rhs();

//Definition of the all matrix of the Darcy System M=(A,B;-B^T,0)
Eigen::MatrixXd A_=A.get_matrix();
Eigen::MatrixXd B_=B.get_matrix();

//Definition of the first row of M (A,B), we call it C
Eigen::MatrixXd C(A_.rows(),A_.cols()+B_.cols());
C<<A_,B_;

//Definition of the second row of M (-B^T,0), we call it E
Eigen::MatrixXd D= Eigen::MatrixXd::Zero(B_T.rows(),B_.cols());
Eigen::MatrixXd E(B_T.rows(),B_T.cols()+D.cols());
E<<-B_T,D;

//Assembling of the matrix M
Eigen::MatrixXd M(C.rows()+E.rows(),C.cols());
M<<C,E;

//Definition of the rhs of the Darcy System. The rhs of the first equation and the one of the second are glued to build the all rhs.
Eigen::VectorXd v1=A.get_rhs();
Eigen::VectorXd v2=B.get_rhs();
Eigen::VectorXd v3(v1.size()+v2.size());
v3<<v1,v2;

//We return the Matrix of the linear system and its rhs.
rhs=v3;
Matrix=M;
}

void Transport_system_esplicit(Eigen::MatrixXd &solution,Eigen::VectorXd &vel,Data_Transport &data)
{
auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,C0]=data;
double h =static_cast<double>(L)/Nx;
double dt=static_cast<double>(T)/Nt;

Eigen::VectorXd IC=Eigen::VectorXd::Zero(Nx);
for (unsigned int i=0;i<Nx/2;++i)
               IC(i)=C0(0.05+i*0.1);
  
Eigen::MatrixXd M(Nx,Nx);
Eigen::VectorXd rhs(Nx);
Eigen::VectorXd sol(Nx);
 
Matrix_F_piu F_p(Nx,Nx);
F_p.assemble_matrix(bc_cond,vel);

Matrix_F_meno F_m(Nx,Nx);
F_m.assemble_matrix(bc_cond,vel);

Matrix_C C(Nx,Nx);
C.assemble_matrix(bc_cond,phi,h,C_in);

M=1/dt*C.get_matrix();

solution.col(0)=IC;

for(unsigned int i=1;i<Nt;i++)
  { 
   rhs=(M-F_p.get_matrix()+F_m.get_matrix())*IC+C.get_rhs();
   solution.col(i)=M.fullPivLu().solve(rhs);
   IC=solution.col(i);
  }
}




void Transport_system_implicit(Eigen::MatrixXd &solution,Eigen::VectorXd &vel,Data_Transport &data)
{
auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,C0]=data;
double h =static_cast<double>(L)/Nx;
double dt=static_cast<double>(T)/Nt;

Eigen::VectorXd IC=Eigen::VectorXd::Zero(Nx);
for (unsigned int i=0;i<Nx/2;++i)
               IC(i)=C0(0.05+i*0.1);
  
Eigen::MatrixXd M(Nx,Nx);
Eigen::VectorXd rhs(Nx);
Eigen::VectorXd sol(Nx);
 
Matrix_F_piu F_p(Nx,Nx);
F_p.assemble_matrix(bc_cond,vel);

Matrix_F_meno F_m(Nx,Nx);
F_m.assemble_matrix(bc_cond,vel);

Matrix_C C(Nx,Nx);
C.assemble_matrix(bc_cond,phi,h,C_in);

M=1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix();

solution.col(0)=IC;

for(unsigned int i=1;i<Nt;i++)
  { 
   rhs=1/dt*C.get_matrix()*IC+C.get_rhs();
   solution.col(i)=M.fullPivLu().solve(rhs);
   IC=solution.col(i);
  }
}
