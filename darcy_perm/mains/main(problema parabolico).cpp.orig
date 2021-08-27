#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>


int main(int argc, char **argv)
{
Data_Darcy data("data.pot");
//std::cout<<"Dati presi"<<std::endl;

Eigen::VectorXi a(5);
Eigen::VectorXd err(5);
a<<100,200,1600,3200,6400;


muparser_fun sol_ex;
sol_ex.set_value("1e6*(1.0-x*x/2)");//"1e6*(1.0-x); caso lineare
//sol_ex.set_value("1e6*(1.0-x)");//"1e6*(1.0-x); caso lineare

Eigen::VectorXd sol;

Eigen::VectorXd exact;

for (unsigned int i=1;i<6;++i)
{
data.Nx=a(i-1);
std::cout<<data.Nx<<std::endl;
double h =static_cast<double>(data.L)/data.Nx;
std::cout<<h<<std::endl;
sol.resize(data.Nx+data.Nx+1);
exact.resize(data.Nx);

//Eigen::MatrixXd M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::SparseMatrix<double> M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::VectorXd rhs(data.Nx+data.Nx+1);
set_Darcy_system(data,M,rhs,h);
//std::cout<<"Sistema Darcy definito"<<std::endl;
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
sol= solver.solve(rhs); 

for(unsigned int j=0;j<data.Nx;++j)
 exact(j)=sol_ex(h/2+j*h);


err(i-1)=(exact-sol.tail(data.Nx)).norm();

std::cout<<"finito "<<err(i-1)<<std::endl;

if(data.Nx==100)
   {Darcy_output_results(sol,data.Nx,data.L);//
    std::cout<<"ECCO "<<data.Nx<<std::endl;
    }

}

std::cout<<std::log(err(3)/err(4))/std::log(2.0)<<std::endl;
std::cout<<std::log(err(2)/err(3))/std::log(2.0)<<std::endl;
std::cout<<std::log(err(1)/err(2))/std::log(2.0)<<std::endl;
std::cout<<std::log(err(0)/err(1))/std::log(2.0)<<std::endl;


Darcy_output_results(sol,data.Nx,data.L);
output_results(exact,data.Nx,data.L);
output_error(err,a);
//std::cout<<"Output salvato"<<std::endl;




return 0;
} 
