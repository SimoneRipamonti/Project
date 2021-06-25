#include "matrix.hpp"
#include "Eigen/Dense"
#include <iostream>
#include <string>

int main()
{
const std::string eccoci="Out";
double c=1.0;
const Eigen::VectorXd vel=-2*Eigen::VectorXd::Ones(5);
Matrix_F_piu Fpiu(4,4);
Fpiu.assemble_matrix(eccoci,c,vel);
Matrix_F_meno Fmeno(4,4);
Fmeno.assemble_matrix(eccoci,c,vel);
std::cout<<vel<<std::endl;
std::cout<<Fpiu.get_matrix()<<std::endl;
std::cout<<" "<<std::endl;
std::cout<<-Fmeno.get_matrix()<<std::endl;
std::cout<<" "<<std::endl;
std::cout<<Fpiu.get_rhs()-Fmeno.get_rhs()<<std::endl;


}
