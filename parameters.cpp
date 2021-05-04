#include "parameters.hpp"
#include "GetPot"
#include <string>

Data::Data(const std::string &filename):K("x"),phi("x+1"),f("x+2")
{
  GetPot file(filename.c_str());
  
  domain_length=file("domain/domain_length",1.);
  const std::string ciao=file("physical_parameters/K_distr","1.0e-10+0*x");
  muparser_fun K(ciao);
  const std::string ecco=file("physical_parameters/phi0","0.2+0*x");
  muparser_fun phi(ecco);
  mu=file("physical_parameters/mu",1.0e-3);
  Q_in=file("BC/Q_in",0.);
  Q_out=file("BC/Q_out",1.0e-1);
  p_in=file("BC/p_in",1.0e6);
  p_out=file("BC/p_out",0.);
  const std::string source=file("Source/f","1.0e-1*x");
  muparser_fun f(source);
  Nx=file("Discretization/Nx",100);
  BC_in=file("BC/in","Pressure");
  BC_out=file("BC/out","Flow");
  dt=file("Time/dt",0.01);
  final_time=file("Time/T",2.0);
  C_in=file("BC_trac/C_in",1.0);
  C_out=file("BC_trac/C_out",1.0);
  bc_cond=file("BC_trac/bc_cond",1.0);
}
