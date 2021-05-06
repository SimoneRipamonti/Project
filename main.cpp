#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "functions.hpp"


int main(int argc, char **argv)
{
Data data("data.pot");
Eigen::MatrixXd M;
Eigen::VectorXd rhs;
set_Darcy_system(data,M,rhs);
Eigen::MatrixXd sol=M.fullPivLu().solve(rhs);
std::cout<<sol<<std::endl;

/*Eigen::VectorXd v(3);
v<<1,2,3;
Eigen::VectorXd a(3);
a<<1,2,3;
std::string ciao="ciao";
output_results(v,a,ciao);
*/



/*Eigen::MatrixXd A=Eigen::MatrixXd::Zero(3,3);
Eigen::MatrixXd B=Eigen::MatrixXd::Zero(3,2);
Eigen::MatrixXd C(A.rows(),A.cols()+B.cols());
C<<A,B;
Eigen::MatrixXd D=Eigen::MatrixXd::Zero(1,5);
Eigen::MatrixXd E(C.rows()+D.rows(),D.cols());
E<<C,D;
std::cout<<E<<std::endl;

Eigen::VectorXd v(3);
v<<1,2,3;
Eigen::VectorXd v2(5);
v2<<4,5,6,7,8;
Eigen::VectorXd v3(v.size()+v2.size());
v3<<v,v2;
std::cout<<v3<<std::endl;
std::cout<<"Hello world!"<<std::endl;
std::cout<<"This is a new added line"<<std::endl;*/

return 0;
} 
