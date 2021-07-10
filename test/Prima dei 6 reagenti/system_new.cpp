//Loop:
//1) Calcolo dei vari rd.
//2) Funzione che mi risolve uno step temporale dei problemi di trasporto e reazione.
//3) Funzione che mi risolve il sistemino non lineare per trovare le concentrazioni effettive. 

#include <Eigen/Dense>
#include <cmath>
#include <string>
#include "functions.hpp"

//Dati
Data_Reaction data_reaction("data.pot");
auto [Ca_0, H_piu_0, HCO3_meno_0, CO2_0, CaSiO3_0, SiO2_0, A, Rate_const, E, R, Temperature, K_eq, K_sol, n, L, Nx, T, Nt, phi]=data_reaction;
std::string bc("In");
double C_in,C_out=0.0;
std::Eigen::VectorXd::Ones(Nx) vel;
double h=static_cast<double>(1/Nx);

//Definizione variabili
Eigen::MatrixXd::Zeros(Nx,Nt) Ca;
Eigen::MatrixXd::Zeros(Nx,Nt) H_piu;
Eigen::MatrixXd::Zeros(Nx,Nt) HCO3_meno;
Eigen::MatrixXd::Zeros(Nx,Nt) CO2;
Eigen::MatrixXd::Zeros(Nx,Nt) CaSiO3;
Eigen::MatrixXd::Zeros(Nx,Nt) SiO2;

Eigen::VectorXd::Zeros(Nx) phi1;
Eigen::VectorXd::Zeros(Nx) phi2;
Eigen::VectorXd::Zeros(Nx) phi3;
Eigen::VectorXd::Zeros(Nx) phi4;
Eigen::VectorXd::Zeros(Nx) phi5;



//Set the initial condition
set_initial_cond(Ca,H_piu,HCO3_meno,CO2,CaSiO3,SiO2,Ca_0,H_piu_0,HCO3_meno_0.,CO2_0,CaSiO3_0,SiO2_0,h);

//Setting for Transport and Reaction Equation
Eigen::MatrixXd M(Nx,Nx);
Eigen::VectorXd rhs(Nx);
assemble_transport(M,rhs,vel,phi,h,Nx);
const auto M_lu=M.fullPivLu();

//Constante per reazione
const double const_r= A*Rate_const*(std::exp(-E/(R*Temperature)));

//Termine di reazione
Eigen::VectorXd::Zeros(Nx) rd;
////////

//Loop Temporale

for(unsigned int step=1; step<T; ++step)
{ 
  compute_phi(phi1,phi2,phi3,phi4,phi5,Ca.col(step-1),H_piu.col(step-1),HCO3_meno.col(step-1),CO2.col(step-1),CaSiO3.col(step-1),SiO2.col(step-1)); //Calcolo le phi
  
  compute_rd(rd,Ca.col(step-1),H_piu.col(step-1),SiO2.col(step-1),const_r,K_sol,n);//Calcolo i termini di reazione
  
  one_step_transport_reaction(phi1,phi2,phi3,phi4,phi5,rd,M_lu,rhs); //Calcolo un passo della Reazione
  
  compute_concentration(Ca,H_piu,HCO3,CO2,CaSiO3,SiO2,phi1,phi2,phi3,phi4,phi5,step,K_eq); //Calcolo le Concentrazioni effettive
}




