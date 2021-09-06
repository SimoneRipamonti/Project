//#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "systems.hpp"
#include <fstream>
#include <exception>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>


//main that simulates only the transport part of the problem
int main(int argc, char **argv)
{

    Data_Transport data_transport("data.pot");//We get the data we need for transport
    Data_linear_decay initial_cond("data.pot");//We get also the linear decay data because inside its definition there is the initial condition of our tracer

    Eigen::VectorXd vel{Eigen::VectorXd::Ones(data_transport.Nx+1)};//transport velocity, we take for semplicity velocity equal to 1[m/s]


    Eigen::MatrixXd Ca(data_transport.Nx,data_transport.Nt+1);//matrix where we store our solution for the tracer (the rows are the spatial steps and the columns are the temporal ones) 

    //We choose between the Esplicit or Implicit scheme
    if(data_transport.method=="Esplicit")
        Transport_system_esplicit(Ca,vel,data_transport,initial_cond);
    else if(data_transport.method=="Implicit")
        Transport_system_implicit(Ca,vel,data_transport,initial_cond);
    else throw std::invalid_argument("Invalid argument: wrong input method, choose  or implicit or esplicit");

    //output results
    output_results_fixed_time("Ca",Ca,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

    output_results_fixed_space("Ca",Ca,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

}
