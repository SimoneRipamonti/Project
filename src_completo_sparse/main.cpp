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

#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>



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
Data_CO2 data_CO2("data.pot");

unsigned int Nx{concentration.get_Nx()};
unsigned int Nt{concentration.get_Nt()};

Eigen::VectorXd vel{0.1*Eigen::VectorXd::Ones(Nx+1)};


Eigen::VectorXd phi1{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi2{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi3{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi4{Eigen::VectorXd::Zero(Nx)};
Eigen::VectorXd phi5{Eigen::VectorXd::Zero(Nx)};



//change_element (Ca);
//Set the initial condition

concentration.set_initial_cond();


//Setting for Transport and Reaction Equation
//Eigen::MatrixXd M(Nx,Nx);
Eigen::SparseMatrix<double> M(Nx,Nx);
Eigen::SparseMatrix<double> rhs(Nx,Nx);
concentration.assemble_transport(M,rhs,vel);



auto [C_in,C_out,bc_cond]=data_CO2;
Eigen::SparseMatrix<double> M_CO2(Nx,Nx);
Eigen::VectorXd rhs_CO2(Nx);


concentration.assemble_transport_CO2(M_CO2,rhs_CO2,vel,C_in,C_out,bc_cond);

Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
// fill A and b;
// Compute the ordering permutation vector from the structural pattern of A
solver.analyzePattern(M); 
// Compute the numerical factorization 
solver.factorize(M); 


Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver2;
// fill A and b;
// Compute the ordering permutation vector from the structural pattern of A
solver2.analyzePattern(M_CO2); 
// Compute the numerical factorization 
solver2.factorize(M_CO2); 

int method{1};

//Termine di reazione
Eigen::VectorXd rd{Eigen::VectorXd::Zero(Nx)};


//Loop Temporale

for(unsigned int i=1; i<Nt; i++)
{ 
  concentration.compute_phi(i-1,phi1,phi2,phi3,phi4,phi5); //Calcolo le phi

  concentration.compute_rd_kp(i-1,rd);//Calcolo i termini di reazione
  
  //concentration.one_step_transport_reaction(phi1,phi2,phi3,phi4,phi5,rd,M,rhs,method,i,solver); //Calcolo un passo della Reazione
 
  concentration.one_step_transport_reaction(phi1,phi2,phi3,phi4,phi5,rd,rhs,rhs_CO2,method,i,solver,solver2);
  
  concentration.compute_concentration(i,phi1,phi2,phi3,phi4,phi5); //Calcolo le Concentrazioni effettive
}


concentration.output_results_fixed_time("Ca");

concentration.output_results_fixed_space("Ca");

concentration.output_results_fixed_space("CaSiO3");

concentration.output_all_reagents(Nx-1);

}
