#include "output.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>



void Darcy_output_results(Eigen::VectorXd &sol,unsigned int Nx,double L) 
{
  // Velocity results to CSV file.
  std::ofstream file1("velocity.csv", std::ofstream::out);
  file1 << "space, velocity" << std::endl;
  Eigen::VectorXd x1=Eigen::VectorXd::LinSpaced(Nx+1,0,L);
  Eigen::VectorXd vel=sol.head(Nx+1);
  
  for (unsigned int i = 0;i<vel.size(); ++i)
    {
      file1<< x1[i] << ", " << vel[i] << std::endl;
    }
  file1.close();
  
  // Pressure results to CSV file.
  std::ofstream file2("pressure.csv", std::ofstream::out);
  file2 << "space, pressure" << std::endl;
  Eigen::VectorXd x2=Eigen::VectorXd::LinSpaced(Nx,0.5,L);
  Eigen::VectorXd pressure=sol.tail(Nx);
  
  for (unsigned int i = 0;i<pressure.size(); ++i)
    {
      file2 << x2[i] << ", " << pressure[i] << std::endl;
    }
  file2.close();
}




void Transport_output_results(Eigen::MatrixXd &value,unsigned int Nx, double L,unsigned int Nt)
{
  std::ofstream file("transport.csv", std::ofstream::out);
  file<< "#space, t0,...,t_Nt-1" << std::endl;
  Eigen::VectorXd x=Eigen::VectorXd::LinSpaced(Nx,0,L);
  for (unsigned int i = 0;i<Nx; ++i)
      {file<< x[i] <<", ";
       for (unsigned int j=0;j<Nt;++j)
            file<<value(i,j)<<", ";
       file<<std::endl;}
  file.close();
}




