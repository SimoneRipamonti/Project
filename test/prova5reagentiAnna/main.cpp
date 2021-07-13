//Loop:
//1) Calcolo dei vari rd.
//2) Funzione che mi risolve uno step temporale dei problemi di trasporto e reazione.
//3) Funzione che mi risolve il sistemino non lineare per trovare le concentrazioni effettive. 

#include <Eigen/Dense>
#include <cmath>
#include <string>
#include "functions.hpp"
#include "muparser_fun.hpp"
#include "output.hpp"
#include "parameters.hpp"
#include <iostream>
#include "matrix.hpp"



int main(){
//Dati
Data_Reaction data_reaction("data.pot");
auto [Ca_0, H_piu_0, HCO3_meno_0, CO2_0, CaSiO3_0, A, Rate_const, E, R, Temperature, K_eq, K_sol, n, L, Nx, T, Nt, phi]=data_reaction;
//std::string bc("In");
//double C_in,C_out=0.0;
Eigen::VectorXd vel{Eigen::VectorXd::Zero(Nx+1)};
double h=static_cast<double>(L/Nx);
double dt=static_cast<double>(T/Nt);



/*std::cout<<"kD1="<<Rate_const<<std::endl;
std::cout<<"Ea="<<E<<std::endl;
std::cout<<"R="<<R<<std::endl;
std::cout<<"Temperature="<<Temperature<<std::endl;
std::cout<<"K1="<<K_eq<<std::endl;
std::cout<<"Ksp1="<<K_sol<<std::endl;
std::cout<<"n="<<n<<std::endl;
std::cout<<"Nx="<<Nx<<std::endl;
std::cout<<"Nt="<<Nt<<std::endl;
std::cout<<"L="<<L<<std::endl;
std::cout<<"T="<<T<<std::endl;
*/
//Definizione variabili
Eigen::MatrixXd Ca{Eigen::MatrixXd::Zero(Nx,Nt)};
Eigen::MatrixXd H_piu{Eigen::MatrixXd::Zero(Nx,Nt)};
Eigen::MatrixXd HCO3_meno{Eigen::MatrixXd::Zero(Nx,Nt)};
Eigen::MatrixXd CO2{Eigen::MatrixXd::Zero(Nx,Nt)};
Eigen::MatrixXd CaSiO3{Eigen::MatrixXd::Zero(Nx,Nt)};
//Eigen::MatrixXd::Zeros(Nx,Nt) SiO2;

Eigen::VectorXd phi1{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi2{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi3{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi4{Eigen::VectorXd::Zero(Nx)};
//Eigen::VectorXd::Zeros(Nx) phi5;



//change_element (Ca);
//Set the initial condition

set_initial_cond(Ca,H_piu,HCO3_meno,CO2,CaSiO3,Ca_0,H_piu_0,HCO3_meno_0,CO2_0,CaSiO3_0,h);


//Setting for Transport and Reaction Equation
Eigen::MatrixXd M(Nx,Nx);
Eigen::MatrixXd rhs(Nx,Nx);
assemble_transport(M,rhs,vel,phi,h,Nx,dt);
//const Eigen::MatrixXd M_lu{M.fullPivLu()};

//Constante per reazione
const double const_r= A*Rate_const*(std::exp(-E/(R*Temperature)));

//Termine di reazione
Eigen::VectorXd rd{Eigen::VectorXd::Zero(Nx)};
////////

//Loop Temporale




for(unsigned int i=1; i<Nt; i++)
{ 
  compute_phi(phi1,phi2,phi3,phi4,Ca.col(i-1),H_piu.col(i-1),HCO3_meno.col(i-1),CO2.col(i-1),CaSiO3.col(i-1)); //Calcolo le phi
  
  compute_rd(rd,Ca.col(i-1),H_piu.col(i-1),const_r,K_sol,n);//Calcolo i termini di reazione
  
  one_step_transport_reaction(phi1,phi2,phi3,phi4,rd,M,rhs); //Calcolo un passo della Reazione
  
  compute_concentration(Ca,H_piu,HCO3_meno,CO2,CaSiO3,phi1,phi2,phi3,phi4,i,K_eq); //Calcolo le Concentrazioni effettive
}



output_results_fixed_time(Ca,Nx,L,Nt);

//output_results_fixed_space("Ca_fixed_space.csv",Ca,Nx,T,Nt);

output_results_fixed_space("CaSiO3_fixed_space.csv",CaSiO3,Nx,T,Nt);

//output_results_fixed_space("HCO3_meno_fixed_space.csv",HCO3_meno,Nx,T,Nt);

//output_results_fixed_space("CO2_fixed_space.csv",CO2,Nx,T,Nt);

//output_results_fixed_space("H_piu_fixed_space.csv",H_piu,Nx,T,Nt);

output_all_reagents(Ca,H_piu,HCO3_meno,CO2,CaSiO3,0,T,Nt);

}
