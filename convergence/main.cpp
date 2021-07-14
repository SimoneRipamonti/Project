#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"

int main(int argc, char **argv)
{
Data_Darcy data("data.pot");
std::cout<<"Dati presi"<<std::endl;
Eigen::MatrixXd M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::VectorXd rhs(data.Nx+data.Nx+1);
set_Darcy_system(data,M,rhs);
std::cout<<"Sistema Darcy definito"<<std::endl;
Eigen::VectorXd sol=M.partialPivLu().solve(rhs);
//Eigen::VectorXd sol=M.fullPivLu().solve(rhs);
std::cout<<"Soluzione Darcy Ottenuta"<<std::endl;


Darcy_output_results(sol,data.Nx,data.L);
std::cout<<"Output salvato"<<std::endl;




return 0;
} 
