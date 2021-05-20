#include "parameters.hpp"
#include "GetPot"
#include <string>

Data_Darcy::Data_Darcy(const std::string &filename)
{
  GetPot file(filename.c_str());
  
  L=file("domain/domain_length",1.);
  const std::string stringa1=file("physical_parameters/K_distr","1.0e-10+0*x");
  K.set_value(stringa1);
  const std::string stringa2=file("physical_parameters/phi0","0.2+0*x");
  phi.set_value(stringa2);
  mu=file("physical_parameters/mu",1.0e-3);
  Q_in=file("BC/Q_in",0.);
  Q_out=file("BC/Q_out",1.0e-1);
  p_in=file("BC/p_in",1.0e6);
  p_out=file("BC/p_out",0.);
  const std::string stringa3=file("Source/f","1.0e-1*x");
  f.set_value(stringa3);
  Nx=file("Discretization/Nx",100);
  BC_in=file("BC/in","Pressure");
  BC_out=file("BC/out","Flow");

}


Data_Transport::Data_Transport(const std::string &filename)
{

  GetPot file(filename.c_str());

  L=file("domain/domain_length",1.);
  const std::string stringa4=file("physical_parameters/phi0","0.2+0*x");
  phi.set_value(stringa4);
  Nx=file("Discretization/Nx",100);
  Nt=file("Time/Nt",100);
  final_time=file("Time/T",2.0);
  C_in=file("BC_trac/C_in",1.0);
  C_out=file("BC_trac/C_out",1.0);
  bc_cond=file("BC_trac/bc_trac","In");
  const std::string stringa5=file("Initial_Condition/IC","2.0*x");
  C0.set_value(stringa5);
}