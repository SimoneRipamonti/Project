#include "parameters.hpp"
#include "GetPot"
#include <string.h>
Data::Data(const std::string &filename)
{
  GetPot file(filename.c_str());
  
  domain_length=file("domain/domain_length",1.);
  K(file(("physical_parameters/K_distr").c_str(),"1.0e-10+0*x"));
  phi(file(("physical_parameters/phi0").c_str(),"0.2+0*x"));
  mu=file("physical_parameters/mu",1.0e-3);
  Q_in=file("BC/Q_in",0.);
  Q_out=file("BC/Q_out",1.0e-1);
  p_in=file("BC/p_in",1.0e6);
  p_out=file("BC/p_out",0.);
  f(file(("Source/f").c_str(),"1.0e-1*x"));
  Nx=file("Discretization/Nx",100);
  BC_in=file("BC/in","Pressure");
  BC_out=file("BC/out","Flow");
}
