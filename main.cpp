//Loop:
//1) Calcolo dei vari rd.
//2) Funzione che mi risolve uno step temporale dei problemi di trasporto e reazione.
//3) Funzione che mi risolve il sistemino non lineare per trovare le concentrazioni effettive. 

#include <Eigen/Dense>
#include <cmath>
#include <string>

#include "muparser_fun.hpp"
#include "parameters.hpp"
#include <iostream>
#include "matrix.hpp"
#include "concentrations.hpp"



int main(){

//Parte su Darcy
/*Data_Darcy data_d("data.pot");

Eigen::MatrixXd M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);
Eigen::VectorXd rhs(data.Nx+data.Nx+1);
set_Darcy_system(data_d,M,rhs);
Eigen::VectorXd sol_darcy=M.fullPivLu().solve(rhs);

Darcy_output_results(sol_darcy,data_d.Nx,data_d.L);
*/



Concentration concentration("data.pot");

unsigned int Nx{concentration.get_Nx()};
unsigned int Nt{concentration.get_Nt()};

Eigen::VectorXd vel{Eigen::VectorXd::Zero(Nx+1)};


Eigen::VectorXd phi1{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi2{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi3{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi4{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi5{Eigen::VectorXd::Zero(Nx)};



//change_element (Ca);
//Set the initial condition

concentration.set_initial_cond();


//Setting for Transport and Reaction Equation
Eigen::MatrixXd M(Nx,Nx);
Eigen::MatrixXd rhs(Nx,Nx);
concentration.assemble_transport(M,rhs,vel);
//const Eigen::MatrixXd M_lu{M.fullPivLu()};

int method{3};


//Termine di reazione
Eigen::VectorXd rd{Eigen::VectorXd::Zero(Nx)};


//Loop Temporale

for(unsigned int i=1; i<Nt; i++)
{ 
  concentration.compute_phi(i,phi1,phi2,phi3,phi4,phi5); //Calcolo le phi
  
  concentration.compute_rd(i,rd);//Calcolo i termini di reazione
  
  concentration.one_step_transport_reaction(phi1,phi2,phi3,phi4,phi5,rd,M,rhs,method,i); //Calcolo un passo della Reazione
  
  concentration.compute_concentration(i,phi1,phi2,phi3,phi4,phi5); //Calcolo le Concentrazioni effettive
}


concentration.output_results_fixed_time("Ca");

concentration.output_results_fixed_space("Ca");

concentration.output_results_fixed_space("CaSiO3");

concentration.output_all_reagents(Nx-1);

}
