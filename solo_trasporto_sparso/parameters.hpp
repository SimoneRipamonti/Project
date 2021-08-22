#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include "muparser_fun.hpp"


//Data for the transport equation
class Data_Transport
{
public:
    explicit Data_Transport(const std::string &filename);

    double L; 
    muparser_fun phi;
    unsigned int Nx;
    unsigned int Nt;//Number of time intervals
    double T;//Time length
    double C_in;//In-section data of the tracer
    double C_out;//Out-section data of the tracer
    std::string bc_cond;//Type of condition
    std::string method;
    muparser_fun C_0;
    
    
    //double lambda;
};


#endif







