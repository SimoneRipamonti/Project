#include <Eigen/Dense>
#include <iostream>
int main()
{
Eigen::VectorXd a(5);
Eigen::VectorXd b(5);
Eigen::VectorXd c(5);

a<<1,
   2,
   3,
   4, 
   5;
b.fill(4);
c.fill(2);

Eigen::VectorXd d(a.size()+b.size()+c.size()+1);
d<<a,b,c,10;
//d=(a.cwiseProduct(b)).cwiseQuotient(c)*3; 

std::cout<<d<<std::endl;



}
