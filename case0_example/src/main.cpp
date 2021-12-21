#include <iostream>
#include <Eigen/Dense>
#include "muparser_fun.hpp"
#include "parameters.hpp"
#include "matrix.hpp"
#include "functions.hpp"
#include <Eigen/LU>
#include "types.hpp"
#include "output.hpp"
#include <cmath>
#include <fstream>
#include <string>

//This main shows the resolution of a Poisson problem, taken as example
int main(int argc, char **argv)
{
    Data_example data("data.pot");
    
    auto &[L,Nx,source]=data;
    
    A_example A(Nx,Nx);
   
    A.assemble_matrix(source,Nx,L);
    
    Solver solver; //Initialization of the solver for the sparse system
    set_solver(A.get_matrix(),solver);
    
    Vector sol{solver.solve(A.get_rhs())};//The Darcy system is solved and the solution is stored in the sol vector

    output_example(sol,Nx,L);//Store the output result in csv files
    
    return 0;
}


