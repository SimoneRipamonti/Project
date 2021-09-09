//#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output_2reag.hpp"
#include "systems_2reag.hpp"
#include <fstream>


int main(int argc, char **argv)
{


    Data_Transport data_transport("data.pot");
    Data_Reaction  data_reaction("data.pot");
    Data_2Reagents data_2reagents("data.pot");

    Eigen::VectorXd vel=0.0*Eigen::VectorXd::Ones(data_transport.Nx+1);
    Eigen::MatrixXd Ca(data_transport.Nx,data_transport.Nt+1);
    Eigen::MatrixXd CaSiO3(data_transport.Nx,data_transport.Nt+1);


    Transport_system_implicit_2_reagents(Ca,CaSiO3,vel,data_transport,data_reaction,data_2reagents);

    output_results_fixed_time_2_reagents(Ca,CaSiO3,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

    output_results_fixed_space_2_reagents(Ca,CaSiO3,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

}
