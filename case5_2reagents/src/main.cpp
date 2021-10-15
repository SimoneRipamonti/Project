//#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output_2reag.hpp"
#include "systems_2reag.hpp"
#include <fstream>

//main that show that simulates the CaSiO3 dissoution reaction
int main(int argc, char **argv)
{

    //We get the data we need from the file
    Data_Transport data_transport("data.pot");
    Data_Reaction  data_reaction("data.pot");
    Data_2Reagents data_2reagents("data.pot");
    
    //For this simulation we take the velocity equal to zero for semplicity, in the next case (case6_all) we will consider also the transport part
    Vector vel{Eigen::VectorXd::Zero(data_transport.Nx+1)};
    
    //Matrices that store the result in space and time for Ca and for CaSiO3 reagents
    Eigen::MatrixXd Ca(data_transport.Nx,data_transport.Nt+1);
    Eigen::MatrixXd CaSiO3(data_transport.Nx,data_transport.Nt+1);

    //We treat implicitly the transport part and explicitly the reactive one
    Transport_system_implicit_2_reagents(Ca,CaSiO3,vel,data_transport,data_reaction,data_2reagents);

    //Output results in csv files
    output_results_fixed_time_2_reagents(Ca,CaSiO3,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

    output_results_fixed_space_2_reagents(Ca,CaSiO3,data_transport.L,data_transport.Nx,data_transport.T,data_transport.Nt);

}
