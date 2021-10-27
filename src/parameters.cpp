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
    Q_in=file("BC_vel/Q_in",0.);
    Q_out=file("BC_vel/Q_out",1.0e-1);
    p_in=file("BC_vel/p_in",1.0e6);
    p_out=file("BC_vel/p_out",0.);
    const std::string stringa3=file("Source/f","1.0e-1*x");
    f.set_value(stringa3);
    Nx=file("Discretization/Nx",100);
    BC_in=file("BC_vel/in","Pressure");
    BC_out=file("BC_vel/out","Flow");

}


Data_Transport::Data_Transport(const std::string &filename)
{

    GetPot file(filename.c_str());

    L=file("domain/domain_length",1.);
    const std::string stringa4=file("physical_parameters/phi0","0.2+0*x");
    phi.set_value(stringa4);
    Nx=file("Discretization/Nx",100);
    Nt=file("Time/Nt",100);
    T=file("Time/T",2.0);
    C_in=file("BC_trac/C_in",1.0);
    C_out=file("BC_trac/C_out",1.0);
    bc_cond=file("BC_trac/bc_trac","In");
    method=file("BC_trac/method","Esplicit");
}

Data_linear_decay::Data_linear_decay(const std::string &filename)
{
    GetPot file(filename.c_str());

    const std::string stringa1=file("BC_trac/C_0","1.0+0.0*x");
    C_0.set_value(stringa1);
    lambda=file("Linear_dec/lambda",1.0);
}



Data_Reaction::Data_Reaction(const std::string &filename)
{
    GetPot file(filename.c_str());

    A=file("Reaction/A",1.);
    Rate_const=file("Reaction/Rate_const",1.);
    E=file("Reaction/E",1.);
    R=file("Reaction/R",1.);
    Temperature=file("Reaction/Temp",1.);
    K_eq=file("Reaction/K_eq",1.);
    K_sol=file("Reaction/K_sol",1.);
    n=file("Reaction/n",1.);
    kp_i=file("Reaction/kp_i",3.98e-13);
    reaction_rate_model=file("Reaction/reaction_rate_model","constant");
    method=file("Reaction/method",1);
}



Data_2Reagents::Data_2Reagents(const std::string &filename)
{


    GetPot file(filename.c_str());

    const std::string stringa1=file("Reaction/Ca_0","1.0+0.0*x");
    Ca_0.set_value(stringa1);
    const std::string stringa5=file("Reaction/CaSiO3_0","1.0+0.0*x");
    CaSiO3_0.set_value(stringa5);
    K_sol=file("Reaction/K_sol",1.);
    ph=file("Reaction/ph",7.);
}



Data_6Reagents::Data_6Reagents(const std::string &filename)
{


    GetPot file(filename.c_str());

    const std::string stringa1=file("Reaction/Ca_0","1.0+0.0*x");
    Ca_0.set_value(stringa1);
    const std::string stringa2=file("Reaction/H_piu_0","1.0+0.0*x");
    H_piu_0.set_value(stringa2);
    const std::string stringa3=file("Reaction/HCO3_meno_0","1.0+0.0*x");
    HCO3_meno_0.set_value(stringa3);
    const std::string stringa4=file("Reaction/CO2_0","1.0+0.0*x");
    CO2_0.set_value(stringa4);
    const std::string stringa5=file("Reaction/CaSiO3_0","1.0+0.0*x");
    CaSiO3_0.set_value(stringa5);
    const std::string stringa6=file("Reaction/SiO2_0","1.0+0.0*x");
    SiO2_0.set_value(stringa6);

}


Data_CO2::Data_CO2(const std::string &filename)
{

    GetPot file(filename.c_str());

    C_in=file("BC_trac/C_in",0.0);
    C_out=file("BC_trac/C_out",0.0);

    bc_cond=file("BC_trac/bc_trac","In");

}


Data_BD::Data_BD(const std::string &filename)
{

    GetPot file(filename.c_str());
    C_in.resize(5,0.);
    C_out.resize(5,0.);
    
    // Ca H+ CO2 SiO2 HCO3-

    C_in[0]=file("BC_trac/Ca/C_in",0.0);
    C_out[0]=file("BC_trac/Ca/C_out",0.0);
    
    C_in[1]=file("BC_trac/Hpiu/C_in",0.0);
    C_out[1]=file("BC_trac/Hpiu/C_out",0.0);
    
    C_in[2]=file("BC_trac/CO2/C_in",0.0);
    C_out[2]=file("BC_trac/CO2/C_out",0.0);
    
    C_in[3]=file("BC_trac/SiO2/C_in",0.0);
    C_out[3]=file("BC_trac/SiO2/C_out",0.0);

    C_in[4]=0;
    C_out[4]=0;

    bc_cond=file("BC_trac/bc_trac","In");
    
   

}

 double Data_BD::getBC(const std::string species,const std::string inout ){
 
 	if (inout=="in")
 		return C_in[s_mapStringValues[species]];
	else
		return C_out[s_mapStringValues[species]];
 		
 }



Data_example::Data_example(const std::string &filename)
{

    GetPot file(filename.c_str());

    L=file("domain/domain_length",1.);
    Nx=file("Discretization/Nx",100);
    
    const std::string stringa=file("Source/f","1.0+0.0*x");
    f.set_value(stringa);
    
}







