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
    //double lambda;
};


class Data_Reaction
{
public:
    explicit Data_Reaction(const std::string &filename);

    double A; //Reaction area
    double Rate_const;//Constant rate
    double E;//Activation energy
    double R;//Gas constant
    double Temperature;//Temperature

};

class Data_2Reagents
{
public:
    explicit Data_2Reagents(const std::string &filename);
    
    muparser_fun Ca_0;
    muparser_fun CaSiO3_0;
    double K_sol;
    double ph;

};



#endif







