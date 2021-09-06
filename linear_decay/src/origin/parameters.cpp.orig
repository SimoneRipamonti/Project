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
    method=file("Transport/method","Esplicit");
    //lambda=file("Reaction/lambda",1.0);
}



Data_linear_decay::Data_linear_decay(const std::string &filename)
{
    GetPot file(filename.c_str());

    const std::string stringa1=file("Transport/C_0","1.0+0.0*x");
    Ca_0.set_value(stringa1);
    
    lambda=file("Reaction/lambda",1.0);
}



