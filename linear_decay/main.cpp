//#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"
#include <fstream>

#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


int main(int argc, char **argv)
{

    Data_Transport data_transport("data.pot");
    Data_linear_decay data_linear_decay("data.pot");

    Eigen::VectorXd vel=1.0*Eigen::VectorXd::Ones(data_transport.Nx+1);
    Eigen::MatrixXd Ca(data_transport.Nx,data_transport.Nt);


    if(data_transport.method=="Esplicit")
        Transport_system_esplicit(Ca,vel,data_transport,data_linear_decay);
    else if(data_transport.method=="Implicit")
        Transport_system_implicit(Ca,vel,data_transport,data_linear_decay);
    else throw std::invalid_argument("Invalid argument: wrong input method, choose  or implicit or esplicit");

    output_results_fixed_space("Ca", Ca, data_transport.Nx, data_transport.T, data_transport.Nt);

    output_results_fixed_time("Ca", Ca, data_transport.Nx, data_transport.L, data_transport.Nt);
}

