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


#endif







