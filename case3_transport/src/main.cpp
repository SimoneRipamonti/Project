//#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output.hpp"
#include "transport_decay.hpp"
#include <fstream>
#include <exception>
#include <Eigen/LU>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

int main(int argc, char **argv)
{

    Data_Transport data_transport("data.pot");//We get the data we need for transport
    Data_linear_decay initial_cond("data.pot");//We get also the linear decay data because inside its definition there is the initial condition of our tracer

    Vector vel{Vector::Ones(data_transport.Nx+1)};//Transport velocity, we take for semplicity velocity equal to 1[m/s]


    Matrix_full Ca(data_transport.Nx,data_transport.Nt+1);//Matrix where we store our solution for the tracer (the rows are the spatial steps and the columns are the temporal ones) 

    //We choose between the Explicit or Implicit scheme
    if(data_transport.method=="Explicit")
        Transport_system_explicit(Ca,vel,data_transport,initial_cond);
    else if(data_transport.method=="Implicit")
        Transport_system_implicit(Ca,vel,data_transport,initial_cond);
    else throw std::invalid_argument("Invalid argument: wrong input method, choose  or implicit or explicit");

    //Output results
    output_results_fixed_time("Ca",Ca,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

    output_results_fixed_space("Ca",Ca,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

}
