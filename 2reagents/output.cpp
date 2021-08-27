#include "output.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>




void output_results_fixed_time_2_reagents(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2, unsigned int Nx, double L,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1("Ca_fixed_time.csv", std::ofstream::out);
    std::ofstream file2("CaSiO3_fixed_time.csv", std::ofstream::out);

    const unsigned int a=Nx/10;
    file1<<"space, ";
    file2<<"space, ";
    for (unsigned int i=0; i<Nt; i+=a)
    {
        file1<<"t="+std::to_string(static_cast<double>(i)/Nt)+", ";
        file2<<"t="+std::to_string(static_cast<double>(i)/Nt)+", ";
    }
    file1<<std::endl;
    file2<<std::endl;

    double h =static_cast<double>(L)/Nx;
    const Eigen::VectorXd x(Eigen::VectorXd::LinSpaced(Nx,h/2,L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";
        file2<< x[i] <<", ";

        for (unsigned int j=0; j<Nt; j+=a)
        {
            file1<<value1(i,j)<<", ";
            file2<<value2(i,j)<<", ";
        }

        file1<<std::endl;
        file2<<std::endl;
    }
    file1.close();
    file2.close();
}



void output_results_fixed_space_2_reagents(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2,unsigned int Nx, double T,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1("Ca_fixed_space.csv", std::ofstream::out);
    std::ofstream file2("CaSiO3_fixed_space.csv", std::ofstream::out);

    const unsigned int a=Nt/10;

    file1<<"time, ";
    file2<<"time, ";
    for (unsigned int i=0; i<Nt; i+=a)
    {
        file1<<"x="+std::to_string(static_cast<double>(i)/Nx)+", ";
        file2<<"x="+std::to_string(static_cast<double>(i)/Nx)+", ";
    }
    file1<<std::endl;
    file2<<std::endl;

    double dt=static_cast<double>(T)/Nt;
    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(Nt,0.0,T-dt));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<Nt; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";
        file2<< t[i]<<", ";

        for (unsigned int j=0; j<Nx; j+=a)
        {
            file1<<value1(j,i)<<", ";
            file2<<value2(j,i)<<", ";
        }
        file1<<std::endl;
        file2<<std::endl;
    }
    file1.close();
    file2.close();
}







