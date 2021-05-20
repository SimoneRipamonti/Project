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
  double final_time;
  double C_in;
  double C_out;
  std::string bc_cond;
  muparser_fun C0;
};
#endif

      


      

      

