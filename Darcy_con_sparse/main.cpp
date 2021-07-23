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
Eigen::VectorXd errp(5);
Eigen::VectorXd errvel(5);
a<<10,20,40,80,160;
//a<<1,2,3,4,10;


muparser_fun sol_ex;
sol_ex.set_value("1e7*(1/(_pi*_pi*16)*sin(_pi*4*x)-0.1*x-1/(4*_pi)*x+0.1)");//"1e6*(1.0-x); caso lineare
//sol_ex.set_value("1e6*(1.0-x)");//"1e6*(1.0-x); caso lineare
//sol_ex.set_value("sin(_pi*2*x)-2*_pi*x");


muparser_fun sol_vel;
//sol_vel.set_value("-2*_pi*10^(-7)*cos(2*_pi*x)+2*_pi*10^(-7)");
//std::cout<<sol_ex(0)<<std::endl;
sol_vel.set_value("-1/(4*_pi)*cos(4*_pi*x)+0.1+1/(4*_pi)");


Eigen::VectorXd sol;

Eigen::VectorXd exact;
Eigen::VectorXd exact_vel;

for (unsigned int i=1;i<6;++i)
{
data.Nx=a(i-1);
std::cout<<data.Nx<<std::endl;
double h =static_cast<double>(data.L)/data.Nx;
std::cout<<h<<std::endl;
sol.resize(data.Nx+data.Nx+1);
exact.resize(data.Nx);
exact_vel.resize(data.Nx+1);

//Eigen::MatrixXd M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::SparseMatrix<double> M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::VectorXd rhs(data.Nx+data.Nx+1);
set_Darcy_system(data,M,rhs,h);

//std::cout<<rhs<<std::endl;
//std::cout<<"fine rhs"<<std::endl;
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
for(unsigned int j=0;j<data.Nx+1;++j)
 exact_vel(j)=sol_vel(j*h); 


errp(i-1)=(exact-sol.tail(data.Nx)).norm();
errvel(i-1)=(exact_vel-sol.head(data.Nx+1)).norm();

std::cout<<"finito "<<errvel(i-1)<<std::endl;

if(data.Nx==80)
   {Darcy_output_results(sol,data.Nx,data.L);//
    std::cout<<"ECCO "<<data.Nx<<std::endl;
    }

}

//std::cout<<std::log(errp(3)/errp(4))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errp(2)/errp(3))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errp(1)/errp(2))/std::log(2.0)<<std::endl;
//std::cout<<std::log(errp(0)/errp(1))/std::log(2.0)<<std::endl;

std::cout<<std::log(errvel(3)/errvel(4))/std::log(2.0)<<std::endl;
std::cout<<std::log(errvel(2)/errvel(3))/std::log(2.0)<<std::endl;
std::cout<<std::log(errvel(1)/errvel(2))/std::log(2.0)<<std::endl;
std::cout<<std::log(errvel(0)/errvel(1))/std::log(2.0)<<std::endl;


Darcy_output_results(sol,data.Nx,data.L);
pressure_exact_result(exact,data.Nx,data.L);
velocity_exact_result(exact_vel,data.Nx,data.L);
output_error(errp,a);
//std::cout<<"Output salvato"<<std::endl;



return 0;
} 
