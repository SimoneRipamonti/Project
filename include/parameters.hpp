#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include "muparser_fun.hpp"

/*!
  *\brief Data for the Darcy system.
  *
  *This class contains variables related to the physical setting of the problem, informations related to the boundary conditions for the Darcy's equations and values related to the spatial discretization
  */
class Data_Darcy
{
public:
    explicit Data_Darcy(const std::string &filename);/*!<Constructor*/

    double L; /*!<Domain length*/
    muparser_fun K;/*!<Soil permeablity (it could change in space)*/
    muparser_fun phi;/*!<Soil porosity (it could change in space)*/
    double mu;/*!<Fluid viscosity*/
    double Q_in;/*!<Inflow velocity*/
    double Q_out;/*!<Outflow velocity*/
    double p_in;/*!<In-section pressure*/
    double p_out;/*!<Out-section pressure*/ 
    muparser_fun f;/*!<External source*/
    unsigned int Nx;/*!<Number of spatial intervals*/
    std::string BC_in;/*!<Type of in-section data*/
    std::string BC_out;/*!Type of out-section data*/
};

/*!
  *\brief Data for the Transport equation.
  *
  *This class contains variables related to the physical setting of the problem, informations related to the boundary conditions for the Transport's equations and values related to the spatial/time discretization 
  */
class Data_Transport
{
public:
    explicit Data_Transport(const std::string &filename);/*!<Constructor*/

    double L;/*!<Domain length*/
    muparser_fun phi;/*!<Soil porosity*/
    unsigned int Nx;/*!<Number of spatial intervals*/
    unsigned int Nt;/*!Number of time intervals*/
    double T;/*!<Time interval*/
    double C_in;/*!<In-section data of the tracer*/
    double C_out;/*!<Out-section data of the tracer*/
    std::string bc_cond;/*!<Type of boundary condition*/
    std::string method;/*!<Type of temporal numerical scheme chosen (Esplicit/Implicit)*/
};


/*!
  *\brief Data for the Linear decay problem.
  *
  *This class contains only information related to the specific linear decay case, so the tracer initial condition and the characteristic linear decay lambda coefficient 
*/
class Data_linear_decay
{
public:
   explicit Data_linear_decay(const std::string &filename);/*!<Constructor*/
  
   muparser_fun C_0;/*!<Tracer initial condition*/
   double lambda;/*!<Linear decay rate*/
};



/*!
  *\brief Data for the reaction setting .
  *
  *This class contains data that describe the physical setting where the reactions take place
*/
class Data_Reaction
{
public:
    explicit Data_Reaction(const std::string &filename);/*!<Constructor*/

    double A; /*!<Surface where the reaction takes place*/
    double Rate_const;/*!<Constant reaction rate*/
    double E;/*!<Activation energy*/
    double R;/*!<Universal gas constant*/
    double Temperature;/*!<Temperature*/
 

};

/*!
  *\brief Data for 2 reagents case.
  *
  *This class contains data that give the information for the following reaction: CaSiO3 + 2H+ <--> Ca2+ + SiO2 + H2O 
*/

class Data_2Reagents
{
public:
    explicit Data_2Reagents(const std::string &filename);/*!<Constructor*/
    
    muparser_fun Ca_0;/*!<Initial condition for Ca*/
    muparser_fun CaSiO3_0;/*!<Initial condition for CaSiO3*/
    double K_sol;/*!<Reaction equilibrium constant*/
    double ph;/*!<Underground water ph*/

};


/*!
  *\brief Data for 6 reagents case.
  *
  *This class contains data that give the information for the following reactions: CaSiO3 + 2H+ <--> Ca2+ + SiO2 + H2O and H2CO3<-->H+ + HCO3-
*/

class Data_6Reagents
{
public:
    explicit Data_6Reagents(const std::string &filename);/*!<Constructor*/
    
    muparser_fun Ca_0;/*!<Initial condition for Ca*/
    muparser_fun H_piu_0;/*!<Initial condition for H_piu*/
    muparser_fun HCO3_meno_0;/*!<Initial condition for HCO3_meno*/
    muparser_fun CO2_0;/*!<Initial condiiton for CO2*/
    muparser_fun CaSiO3_0;/*!<Initial condition for the Wollastonite*/
    muparser_fun SiO2_0;/*!<Initial condition for SiO2*/

    double K_eq;/*!<Equilibrium constant for the faster reaction H2CO3<-->H+ + HCO3-*/
    double K_sol;/*!<Equilibrium constant for the dissolution reaction  CaSiO3 + 2H+ <--> Ca2+ + SiO2 + H2O*/
    double n;/*!<Exponent for [H]+ in the reaction rate computation*/

    double kp_i;/*!<Precipitation rate constant for the reaction rate computation*/

    int method;/*!<Numerical scheme adopted for the reaction term (Explicit, Heun, Predictor-Corrector*/

};


/*!
  *\brief Data for CO2 inflow
  *
  *This class is specific for the particular test case analyzed: we need specific information for the  CO2 constant inflow
*/
class Data_CO2
{

public:
     explicit Data_CO2(const std::string &filename);/*!<Constructor*/

     double C_in;/*!<CO2 inflow condition*/
     double C_out;/*!<CO2 outflow condition*/
     std::string bc_cond;/*!<Type of boundary condition (Inflow or outflow bc)*/

};



#endif







