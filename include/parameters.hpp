#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include "muparser_fun.hpp"

//Data for the Darcy System
class Data_Darcy
{
public:
    explicit Data_Darcy(const std::string &filename);

    double L; //Domain length
    muparser_fun K;//Permeablity
    muparser_fun phi;//Porosity
    double mu;//Viscosity
    double Q_in;//Inflow velocity
    double Q_out;//Outflow velocity
    double p_in;//In-section pressure data
    double p_out;//Out-section pressure data
    muparser_fun f;//External source
    unsigned int Nx;//Number of spatial intervals
    std::string BC_in;//Type of in-section data
    std::string BC_out;//Type of out-section data
};


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
    //double lambda;
};

class Data_linear_decay
{
public:
   explicit Data_linear_decay(const std::string &filename);
  
   muparser_fun C_0;
   double lambda;//Linear decay rate
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


class Data_6Reagents
{
public:
    explicit Data_6Reagents(const std::string &filename);
    
    //Initial reagents data
    muparser_fun Ca_0;
    muparser_fun H_piu_0;
    muparser_fun HCO3_meno_0;
    muparser_fun CO2_0;
    muparser_fun CaSiO3_0;
    muparser_fun SiO2_0;

    double K_eq;//Equilibrium constant for the faster reaction H2CO3<-->H+ + HCO3-
    double K_sol;//Equilibrium constant for the dissolution reaction  CaSiO3 + 2H+ <--> Ca2+ + SiO2 + H2O
    double n;

    double kp_i;//Precipitation rate constant

    int method;

};


//Since for the complex case we have a costant CO2 inflow
class Data_CO2
{

public:
     explicit Data_CO2(const std::string &filename);

     double C_in;
     double C_out;
     std::string bc_cond;

};



#endif







