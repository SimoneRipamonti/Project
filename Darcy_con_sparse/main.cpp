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

Eigen::VectorXi a(5); //vector that contains the numbers of steps
Eigen::VectorXd errp(5);//vector pressure error
Eigen::VectorXd errvel(5);//vector velocity error

a<<10,20,40,80,160;//Spatial steps


muparser_fun sol_ex; //muparser_fun that contains the exact solution for pressure
sol_ex.set_value("1e7*(1/(_pi*_pi*16)*sin(_pi*4*x)-0.1*x-1/(4*_pi)*x+0.1)");//"1e6*(1.0-x); caso lineare
//sol_ex.set_value("1e6*(1.0-x)");//"1e6*(1.0-x); caso lineare
//sol_ex.set_value("sin(_pi*2*x)-2*_pi*x");


muparser_fun sol_vel; //muparser_fun that contains the exact solution for velocity
//sol_vel.set_value("-2*_pi*10^(-7)*cos(2*_pi*x)+2*_pi*10^(-7)");
//std::cout<<sol_ex(0)<<std::endl;
sol_vel.set_value("-1/(4*_pi)*cos(4*_pi*x)+0.1+1/(4*_pi)");


Eigen::VectorXd sol;//vector that stores the numerical solution

Eigen::VectorXd exact; //vector that stores the exact solution for pressure
Eigen::VectorXd exact_vel; //vector that stores the exact solution for velocity

//Loop where at each iterate it is changed the number of steps
for (unsigned int i=1;i<6;++i) 
{
data.Nx=a(i-1);

double h =static_cast<double>(data.L)/data.Nx; //space step
std::cout<<"h="<<h<<std::endl;

sol.resize(data.Nx+data.Nx+1); //solution vectors are resized
exact.resize(data.Nx);
exact_vel.resize(data.Nx+1);

Eigen::SparseMatrix<double> M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);//Initialization of the big matrix for the Darcy system
Eigen::VectorXd rhs(data.Nx+data.Nx+1);//Initialization of the rhs of Darcy
set_Darcy_system(data,M,rhs,h);//Definition of the Darcy system Mx=rhs


Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver; //Initialization of the Solver for the sparse system

solver.analyzePattern(M); // Compute the ordering permutation vector from the structural pattern of A

solver.factorize(M); // Compute the numerical factorization 

sol= solver.solve(rhs);//The Darcy system is solved and the solution is stored in the sol vector  

////The exact solution vectors are filled
for(unsigned int j=0;j<data.Nx;++j)
 exact(j)=sol_ex(h/2+j*h);
for(unsigned int j=0;j<data.Nx+1;++j)
 exact_vel(j)=sol_vel(j*h); 

//Errors computation
errp(i-1)=(exact-sol.tail(data.Nx)).norm();
errvel(i-1)=(exact_vel-sol.head(data.Nx+1)).norm();

}

std::cout<<"Order of convergence for pressure:"<<std::log(errp(3)/errp(4))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errp(2)/errp(3))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errp(1)/errp(2))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errp(0)/errp(1))/std::log(2.0)<<std::endl;


std::cout<<"Order of convergence for velocity:"<<std::log(errvel(3)/errvel(4))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errvel(2)/errvel(3))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errvel(1)/errvel(2))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errvel(0)/errvel(1))/std::log(2.0)<<std::endl;


//Plot of the output results
Darcy_output_results(sol,data.Nx,data.L);
pressure_exact_result(exact,data.Nx,data.L);
velocity_exact_result(exact_vel,data.Nx,data.L);
output_error(errp,a);




return 0;
} 
