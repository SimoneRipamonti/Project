#include "functions.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>

void set_Darcy_system(Data &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs)
{
auto &[domain_length, K, phi, mu, Q_in, Q_out, p_in, p_out, f, Nx, BC_in, BC_out, dt, final_time, C_in, C_out, bc_cond]=data;

double h =static_cast<double>(domain_length)/Nx;


Matrix_A A(Nx+1,Nx+1,K,h,mu,BC_in,BC_out,p_in,Q_out);
A.set_matrix();
A.set_BC();
A.set_rhs();

Matrix_B B(Nx+1,Nx,BC_in,BC_out,f,h);//è diversa dalle note, la definisco come la B^T delle note
B.set_matrix();
Eigen::MatrixXd B_T=(B.get_matrix()).transpose();
B.set_BC();
B.set_rhs();

Eigen::MatrixXd A_=A.get_matrix();
Eigen::MatrixXd B_=B.get_matrix();


Eigen::MatrixXd C(A_.rows(),A_.cols()+B_.cols());
C<<A_,B_;

Eigen::MatrixXd D= Eigen::MatrixXd::Zero(B_T.rows(),B_.cols());
Eigen::MatrixXd E(B_T.rows(),B_T.cols()+D.cols());
E<<-B_T,D;

Eigen::MatrixXd M(C.rows()+E.rows(),C.cols());
M<<C,E;

Eigen::VectorXd v1=A.get_rhs();
Eigen::VectorXd v2=B.get_rhs();
Eigen::VectorXd v3(v1.size()+v2.size());
v3<<v1,v2;

rhs=v3;
Matrix=M;
}

void output_results(Eigen::VectorXd &sol,unsigned int Nx,double L) 
{
  // Output results to CSV file.
  std::ofstream file1("velocity.csv", std::ofstream::out);
  file1 << "space, velocity" << std::endl;
  Eigen::VectorXd x1=Eigen::VectorXd::LinSpaced(Nx+1,0,L);
  Eigen::VectorXd vel=sol.head(Nx+1);
  
  for (unsigned int i = 0;i<vel.size(); ++i)
    {
      file1<< x1[i] << ", " << vel[i] << std::endl;
    }
  file1.close();
  
  std::ofstream file2("pressure.csv", std::ofstream::out);
  file2 << "space, pressure" << std::endl;
  Eigen::VectorXd x2=Eigen::VectorXd::LinSpaced(Nx,0.5,L);
  Eigen::VectorXd pressure=sol.tail(Nx);
  
  for (unsigned int i = 0;i<pressure.size(); ++i)
    {
      file2 << x2[i] << ", " << pressure[i] << std::endl;
    }
  file2.close();

 
  // Plot results.
  /*Gnuplot gp;
  gp << "set xlabel 'Space'; set ylabel 'value' ; set key center "
        "right; plot "
     << gp.file1d(std::tie(space, value))
     << "with line linewidth 2 title 'Value'" << std::endl;*/
}





