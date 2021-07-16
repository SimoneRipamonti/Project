#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"
#include <Eigen/LU>
#include<Eigen/IterativeLinearSolvers>
#include<Eigen/Sparse>
#include<Eigen/SparseLU>


int main(int argc, char **argv)
{
Data_Darcy data("data.pot");
std::cout<<"Dati presi"<<std::endl;

double h =static_cast<double>(data.L)/data.Nx;

Eigen::VectorXd exact(data.Nx);


muparser_fun sol_ex;
//sol_ex.set_value("1e6*(1.0-x*x/2)");//"1e6*(1.0-x); caso lineare
sol_ex.set_value("1e6*(1.0-x)");//"1e6*(1.0-x); caso lineare

for(unsigned int j=0;j<data.Nx;++j)
 exact(j)=sol_ex(h/2+j*h);
std::cout<<exact<<std::endl;



//Eigen::MatrixXd M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::SparseMatrix<double> M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::VectorXd rhs(data.Nx+data.Nx+1);
set_Darcy_system(data,M,rhs,h);
std::cout<<"Sistema Darcy definito"<<std::endl;
//Eigen::VectorXd sol=M.partialPivLu().solve(rhs);
//Eigen::VectorXd sol=M.fullPivLu().solve(rhs);
//Eigen::VectorXd sol=M.lu().solve(rhs);

Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
// fill A and b;
// Compute the ordering permutation vector from the structural pattern of A
solver.analyzePattern(M); 
// Compute the numerical factorization 
solver.factorize(M); 
//Use the factors to solve the linear system 
Eigen::VectorXd sol= solver.solve(rhs); 

std::cout<<sol.tail(data.Nx)<<std::endl;

std::cout<<"Soluzione Darcy Ottenuta"<<std::endl;




Darcy_output_results(sol,data.Nx,data.L);
std::cout<<"Output salvato"<<std::endl;






return 0;
} 
