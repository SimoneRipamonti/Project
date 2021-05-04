#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <string>
#include "muparser_fun.hpp"
class Data
{
public:
      explicit Data(const std::string &filename);
      //Data(const Data &data);//da implementare

      double domain_length;
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
      double dt;
      double final_time;
      double C_in;
      double C_out;
      std::string bc_cond;
};
#endif

      


      

      

