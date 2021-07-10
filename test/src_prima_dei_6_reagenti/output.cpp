#include "output.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>



void Darcy_output_results(Eigen::VectorXd &sol,unsigned int Nx,double L) //The solution vector is built by the vector velocity (for the first Nx+1 values) and the vector pressure (for the last Nx values)
{
    // Velocity results to CSV file.
    std::ofstream file1("velocity.csv", std::ofstream::out);
    file1 << "space, velocity" << std::endl;
    Eigen::VectorXd x1=Eigen::VectorXd::LinSpaced(Nx+1,0,L);//Definition of the space vector
    Eigen::VectorXd vel=sol.head(Nx+1);//Definition of the vector velocity solution

    for (unsigned int i = 0; i<vel.size(); ++i) //Loop to save the velocity values with respect to the spatial position on the CSV file
        file1<< x1[i] << ", " << vel[i] << std::endl;

    file1.close();

    // Pressure results to CSV file.
    std::ofstream file2("pressure.csv", std::ofstream::out);
    file2 << "space, pressure" << std::endl;
    Eigen::VectorXd x2=Eigen::VectorXd::LinSpaced(Nx,0.05,L-0.05);//Definition of the space vecotr (Pressure values are stored in the middle of the cell)
    Eigen::VectorXd pressure=sol.tail(Nx);//Definition of the pressure vector solution

    for (unsigned int i = 0; i<pressure.size(); ++i) //Loop to save the pressure values with respect to the spatial position on the CSV file
        file2 << x2[i] << ", " << pressure[i] << std::endl;

    file2.close();
}


void Transport_output_results_fixed_time(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2, unsigned int Nx, double L,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1("Ca_fixed_time.csv", std::ofstream::out);
    std::ofstream file2("CaSiO3_fixed_time.csv", std::ofstream::out);
    file1<< "#space, t0,...,t_Nt-1" << std::endl;
    file2<< "#space, t0,...,t_Nt-1" << std::endl;
    
    double h =static_cast<double>(L)/Nx;
    const Eigen::VectorXd x(Eigen::VectorXd::LinSpaced(Nx,h/2,L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)
    
    for (unsigned int i = 0; i<Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";
        file2<< x[i] <<", ";
        
        for (unsigned int j=0; j<Nt; ++j)
            {file1<<value1(i,j)<<", ";
             file2<<value2(i,j)<<", ";  }

        file1<<std::endl;
        file2<<std::endl;
    }
    file1.close();
    file2.close();
}


void Transport_output_results_fixed_space(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2,unsigned int Nx, double T,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1("Ca_fixed_space.csv", std::ofstream::out);
    std::ofstream file2("CaSiO3_fixed_space.csv", std::ofstream::out);
    file1<< "#time, x0,...,x_Nx-1" << std::endl;
    file2<< "#time, x0,...,x_Nx-1" << std::endl;
    
    double dt=static_cast<double>(T)/Nt;
    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(Nt,0.0,T-dt));//Definition of the space vector (Concnetration values are stored in the middle of the cell)
    
    for (unsigned int i = 0; i<Nt; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";
        file2<< t[i]<<", ";
  
        for (unsigned int j=0; j<Nx; ++j)
            {file1<<value1(j,i)<<", ";
             file2<<value2(j,i)<<", ";
            }
        file1<<std::endl;
        file2<<std::endl;
    }
    file1.close();
    file2.close();
}




