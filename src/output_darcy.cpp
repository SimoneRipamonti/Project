#include "output_darcy.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>
#include <vector>


void Darcy_output_results(Vector &sol,unsigned int Nx,double L) //The solution vector is built by the vector velocity (for the first Nx+1 values) and the vector pressure (for the last Nx values)
{
    // Velocity results to CSV file.
    std::ofstream file1("velocity"+std::to_string(Nx)+".csv", std::ofstream::out);
    file1 << "space, Nx="+std::to_string(Nx)<< std::endl;
    Vector x1{Vector::LinSpaced(Nx+1,0,L)};//Definition of the space vector
    Vector vel{sol.head(Nx+1)};//Definition of the vector velocity solution

    for (unsigned int i = 0; i<vel.size(); ++i) //Loop to save the velocity values with respect to the spatial position on the CSV file
        file1<< x1[i] << ", " << vel[i] << std::endl;

    file1.close();

    const double h =static_cast<double>(L)/Nx;

    // Pressure results to CSV file.
    std::ofstream file2("pressure"+std::to_string(Nx)+".csv", std::ofstream::out);
    file2 << "space, Nx="+std::to_string(Nx)<< std::endl;
    Vector x2{Vector::LinSpaced(Nx,h,L-h)};//Definition of the space vecotr (Pressure values are stored in the middle of the cell)
    Vector pressure{sol.tail(Nx)};//Definition of the pressure vector solution

    for (unsigned int i = 0; i<pressure.size(); ++i) //Loop to save the pressure values with respect to the spatial position on the CSV file
        file2 << x2[i] << ", " << pressure[i] << std::endl;

    file2.close();
}


void pressure_exact_result(Vector &value1, unsigned int Nx, double L)
{
    //Concentration value results to CSV file.
    std::ofstream file1("exact_pressure.csv", std::ofstream::out);
    file1<< "space, exact pressure" << std::endl;

    const double h =static_cast<double>(L)/Nx;
    const Vector x(Vector::LinSpaced(Nx,h/2,L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";


        file1<<value1(i)<<", ";


        file1<<std::endl;

    }
    file1.close();

}


void velocity_exact_result(Vector &value1, unsigned int Nx, double L)
{
    //Concentration value results to CSV file.
    std::ofstream file1("exact_velocity.csv", std::ofstream::out);
    file1<< "space, exact velocity " << std::endl;

    const Vector x(Vector::LinSpaced(Nx+1,0,L));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<Nx+1; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";


        file1<<value1(i)<<", ";


        file1<<std::endl;

    }
    file1.close();

}


void output_error(Vector &value1, Eigen::VectorXi &Nx, const std::string &name)
{
    //Concentration value results to CSV file.
    std::ofstream file1("error_"+name+".csv", std::ofstream::out);
    file1<< "Nx, Error "+name+", p = 1, p = 1.5" << std::endl;


    for (unsigned int i = 0; i<value1.size(); ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< Nx(i) <<", ";


        file1<<value1(i)<<", ";
        
        file1<<1./Nx(i)<<", ";

        file1<<1./(std::pow(Nx(i),1.5))<<", ";
      
        file1<<std::endl;

    }
    file1.close();

}







