#include "output.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>





void output_results_fixed_time(Eigen::MatrixXd &value1, unsigned int Nx, double L,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1("Ca_fixed_time.csv", std::ofstream::out);
 
    file1<< "#space, t0,...,t_Nt-1" << std::endl;
    
    double h =static_cast<double>(L)/Nx;
    const Eigen::VectorXd x(Eigen::VectorXd::LinSpaced(Nx,h/2,L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)
    
    for (unsigned int i = 0; i<Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";
     
        for (unsigned int j=0; j<Nt; ++j)
            {file1<<value1(i,j)<<", ";
            }

        file1<<std::endl;
        
    }
    file1.close();

}


void output_results_fixed_space(const std::string& title, const Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1(title, std::ofstream::out);
    file1<< "time, x0,...,x_Nx-1" << std::endl;
   
    
    double dt=static_cast<double>(T)/Nt;
    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(Nt,0.0,T-dt));//Definition of the space vector (Concnetration values are stored in the middle of the cell)
    
    for (unsigned int i = 0; i<Nt; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";
       
        for (unsigned int j=0; j<Nx; ++j)
            {file1<<value1(j,i)<<", ";
       
            }
        file1<<std::endl;
      
    }
    file1.close();

}

void output_all_reagents(const Eigen::MatrixXd &Ca,const Eigen::MatrixXd &H_piu,const Eigen::MatrixXd &CaSiO3,const Eigen::MatrixXd &CO2,const Eigen::MatrixXd &HCO3_meno, unsigned int j, double T, unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1("all.csv", std::ofstream::out);
    file1<< "time, Ca, H_piu, HCO3_meno, CO2, CaSiO3" << std::endl;
   
    
    double dt=static_cast<double>(T)/Nt;
    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(Nt,0.0,T-dt));//Definition of the space vector (Concnetration values are stored in the middle of the cell)
    
    for (unsigned int i = 0; i<Nt; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";
       
        file1<<Ca(j,i)<<", ";
        file1<<H_piu(j,i)<<", ";
        file1<<CaSiO3(j,i)<<", ";
        file1<<CO2(j,i)<<", ";
        file1<<HCO3_meno(j,i)<<", ";
        
        file1<<std::endl;
      
    }
    file1.close();

}





/*
void Transport_output_results_fixed_space(Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt)
{
    //Concentration value results to CSV file.
    const std::string title="Ca_fixed_space.csv";
    std::ofstream file1(title, std::ofstream::out);
    file1<< "time, x0,...,x_Nx-1" << std::endl;
   
    
    double dt=static_cast<double>(T)/Nt;
    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(Nt,0.0,T-dt));//Definition of the space vector (Concnetration values are stored in the middle of the cell)
    
    for (unsigned int i = 0; i<Nt; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";
       
        for (unsigned int j=0; j<Nx; ++j)
            {file1<<value1(j,i)<<", ";
       
            }
        file1<<std::endl;
      
    }
    file1.close();

}

*/





