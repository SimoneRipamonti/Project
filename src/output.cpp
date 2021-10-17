#include "output.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>



void output_results_fixed_time(const std::string& title, const Matrix_full &value1, double L, unsigned int Nx, double T,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1( title+"_fixed_time.csv", std::ofstream::out);
    const double dt=static_cast<double>(T)/Nt;
    const double h =static_cast<double>(L)/Nx;

    const unsigned int a=Nt/10; //we print the result on the file each "a" seconds
    file1<<"space, ";
    for (unsigned int i=0.; i<Nt+1; i+=a)
        file1<<"t="+std::to_string(i*dt)+"s, ";
    file1<<std::endl;

    
    const Vector x(Vector::LinSpaced(Nx,h/2,L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)
  
    for (unsigned int i = 0; i<Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";

        for (unsigned int j=0; j<Nt+1; j+=a)
        {
            file1<<value1(i,j)<<", ";
        }

        file1<<std::endl;

    }
    file1.close();

}


void output_results_fixed_space(const std::string& title, const Matrix_full &value1, double L, unsigned int Nx, double T,unsigned int Nt)
{
    //Concentration value results to CSV file.
    std::ofstream file1(title+"_fixed_space.csv", std::ofstream::out);
    const double h =static_cast<double>(L)/Nx;

    const unsigned int a=Nx/10;//we print the result on the file each "a" meters

    file1<<"time, ";
    for (unsigned int i=0; i<Nx; i+=a)
        file1<<"x="+std::to_string(h/2+i*h)+"m, ";
    file1<<"x="+std::to_string(L-h/2)+"m, ";
    file1<<std::endl;

    
    const Vector t(Vector::LinSpaced(Nt+1,0.0,T));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<Nt+1; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";

        for (unsigned int j=0; j<Nx; j+=a)
        {
            file1<<value1(j,i)<<", ";

        }
        file1<<value1(Nx-1,i)<<", ";
        file1<<std::endl;

    }
    file1.close();

}


