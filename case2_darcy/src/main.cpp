#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output_darcy.hpp"
#include "darcy.hpp"
#include "functions.hpp"
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>


//This main shows the behaviour of the Darcy code with a permeability that is not continous 
int main(int argc, char **argv)
{
    Data_Darcy data("data.pot");//We get the data from a file

    double h =static_cast<double>(data.L)/data.Nx; //space step

    Vector sol(data.Nx+data.Nx+1); //solution vectors are resized

    Matrix M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);//Initialization of the big matrix for the Darcy system
    Vector rhs(data.Nx+data.Nx+1);//Initialization of the rhs of Darcy
    set_Darcy_system(data,M,rhs,h);//Definition of the Darcy system Mx=rhs
    
    Solver solver; //Initialization of the solver for the sparse system
    set_solver(M,solver);
    
    sol= solver.solve(rhs);//The Darcy system is solved and the solution is stored in the sol vector
    
    std::cout<<"M="<<std::endl;
    std::cout<<M<<std::endl;
 
    Darcy_output_results(sol,data.Nx,data.L);//Store the output result in csv files





    return 0;
}
