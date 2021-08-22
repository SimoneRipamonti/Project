#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include "muparser_fun.hpp"

//Data for the Darcy System
class Data_Darcy
{
public:
    explicit Data_Darcy(const std::string &filename);

    double L;
    muparser_fun K;
    muparser_fun phi;
    double mu;
    double Q_in;
    double Q_out;
    double p_in;
    double p_out;
    muparser_fun f;
    unsigned int Nx;
    std::string BC_in;
    std::string BC_out;
};


//Data for the transport equation
class Data_Transport
{
public:
    explicit Data_Transport(const std::string &filename);

    double L;
    muparser_fun phi;
    unsigned int Nx;
    unsigned int Nt;
    double T;
    double C_in;
    double C_out;
    std::string bc_cond;
    //double lambda;
};

class Data_6Reagents
{
public:
    explicit Data_6Reagents(const std::string &filename);
    
    muparser_fun Ca_0;
    muparser_fun H_piu_0;
    muparser_fun HCO3_meno_0;
    muparser_fun CO2_0;
    muparser_fun CaSiO3_0;
    muparser_fun SiO2_0;
    double K_eq;
    double K_sol;
    double n;

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



class Data_Reaction
{
public:
    explicit Data_Reaction(const std::string &filename);

    double A;
    double Rate_const;
    double E;
    double R;
    double Temperature;


};

#endif







