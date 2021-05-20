#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"
#include <fstream>


int main(int argc, char **argv)
{
Data_Transport data("data.pot");

Eigen::VectorXd vel=1*Eigen::VectorXd::Ones(data.Nx+1);
Eigen::MatrixXd solution(data.Nx,data.Nt);

Transport_system_esplicit(solution,vel,data);

Transport_output_results(solution,data.Nx,data.L,data.Nt);
}
