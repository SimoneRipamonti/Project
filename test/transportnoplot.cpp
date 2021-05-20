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

Eigen::VectorXd vel=1*Eigen::VectorXd::Ones(Nx+1);
double h =static_cast<double>(domain_length)/Nx;
unsigned int Nt=final_time/dt;
std::cout<<Nt<<std::endl;
std::cout<<Nx<<std::endl;
Eigen::VectorXd C0=Eigen::VectorXd::Zero(Nx);
//Eigen::MatrixXd solution(Nx,Nt);
//solution.col(0)=C0;
for (unsigned int i=0;i<Nx/2;++i)
               C0(i)=1.;
  
Eigen::MatrixXd M(Nx,Nx);
Eigen::VectorXd rhs(Nx);
Eigen::VectorXd sol(Nx);
 
Matrix_F_piu F_p(Nx,Nx,bc_cond,vel);
Matrix_F_meno F_m(Nx,Nx,bc_cond,vel);
F_p.set_matrix();
F_p.set_rhs();
F_p.set_BC();
F_m.set_matrix();
F_m.set_rhs();
F_m.set_BC();
Matrix_C C(Nx,Nx,bc_cond,phi,h,C_in);
C.set_matrix();
C.set_rhs();
C.set_BC();
M=1/dt*C.get_matrix();

std::cout<<-F_p.get_matrix()+F_m.get_matrix()<<std::endl;

//std::cout<<M<<std::endl;
std::cout<<C0<<std::endl;
std::cout<<" "<<std::endl;


std::ofstream file2("output.csv", std::ofstream::out);
file2 << "space,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9" << std::endl;
Eigen::VectorXd x2=Eigen::VectorXd::LinSpaced(Nx,0.5,L);

for (unsigned int i = 0;i<C0.size(); ++i)
    file2 << x2[i] << ", " << C0[i] <<","<< std::endl;

for(unsigned int i=1;i<Nt;i++)
{ 
//rhs=1/dt*C.get_matrix()*C0+C.get_rhs()+F_m.get_rhs()+F_p.get_rhs();
rhs=(1/dt*C.get_matrix()-F_p.get_matrix()+F_m.get_matrix())*C0;
sol=M.fullPivLu().solve(rhs);
std::cout<<' '<<std::endl;
std::cout<<sol<<std::endl;
C0=sol;
}





return 0;
} 
