#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output_darcy.hpp"
#include "darcy.hpp"
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>


int main(int argc, char **argv)
{
    Data_Darcy data("data.pot");

    double h =static_cast<double>(data.L)/data.Nx; //space step

    Eigen::VectorXd sol(data.Nx+data.Nx+1); //solution vectors are resized

    Eigen::SparseMatrix<double> M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);//Initialization of the big matrix for the Darcy system
    Eigen::VectorXd rhs(data.Nx+data.Nx+1);//Initialization of the rhs of Darcy
    set_Darcy_system(data,M,rhs,h);//Definition of the Darcy system Mx=rhs


    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver; //Initialization of the Solver for the sparse system

    solver.analyzePattern(M); // Compute the ordering permutation vector from the structural pattern of A

    solver.factorize(M); // Compute the numerical factorization

    sol= solver.solve(rhs);//The Darcy system is solved and the solution is stored in the sol vector

 
    Darcy_output_results(sol,data.Nx,data.L);





    return 0;
}
