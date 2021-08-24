#include "parameters.hpp"
#include "GetPot"
#include <string>



Data_Transport::Data_Transport(const std::string &filename)
{

    GetPot file(filename.c_str());

    L=file("domain/domain_length",1.);
    const std::string stringa4=file("physical_parameters/phi0","0.2+0*x");
    phi.set_value(stringa4);
    Nx=file("Discretization/Nx",100);
    Nt=file("Time/Nt",100);
    T=file("Time/T",2.0);
    C_in=file("Transport/C_in",1.0);
    C_out=file("Transport/C_out",1.0);
    bc_cond=file("Transport/bc_trac","In");

}

Data_Reaction::Data_Reaction(const std::string &filename)
{
    GetPot file(filename.c_str());

    A=file("Reaction/A",1.);
    Rate_const=file("Reaction/Rate_const",1.);
    E=file("Reaction/E",1.);
    R=file("Reaction/R",1.);
    Temperature=file("Reaction/Temp",1.);
}


Data_2Reagents::Data_2Reagents(const std::string &filename)
{


    GetPot file(filename.c_str());

    const std::string stringa1=file("Reaction/Ca_0","1.0+0.0*x");
    Ca_0.set_value(stringa1);
    const std::string stringa5=file("Reaction/CaSiO3_0","1.0+0.0*x");
    CaSiO3_0.set_value(stringa5);
    K_sol=file("Reaction/K_sol",1.);
    ph=file("Reaction/ph",7.);
}






