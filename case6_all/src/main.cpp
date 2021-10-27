//Loop:
//1) Calcolo dei vari rd.
//2) Funzione che mi risolve uno step temporale dei problemi di trasporto e reazione.
//3) Funzione che mi risolve il sistemino non lineare per trovare le concentrazioni effettive.

#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <iostream>

#include "darcy.hpp"
#include "muparser_fun.hpp"
#include "parameters.hpp"
#include "matrix.hpp"
#include "concentrations.hpp"
#include "functions.hpp"

#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>



int main()
{
    //Darcy part for computing velocity
    Data_Darcy data("data.pot");

    Vector vel(data.Nx+1);

    Darcy_velocity(data,vel);
    
    //vel=Eigen::VectorXd::Zero(data.Nx+1);

    Concentration concentration("data.pot");
    
    //Set the initial condition for the reagents
    concentration.set_initial_cond();

    unsigned int Nx{concentration.get_Nx()};
    unsigned int Nt{concentration.get_Nt()};
    
    //Eigen::VectorXd vel{6.67e-9*Eigen::VectorXd::Ones(Nx+1)};

    Solver  solver, solver1;
    Matrix M_rhs(Nx,Nx);
    Vector rhs_psi1(Nx);
    Vector rhs_psi2(Nx);
    Vector rhs_psi3(Nx);
    Vector rhs_psi5(Nx);
    concentration.define_transport_solver(solver, solver1, M_rhs, rhs_psi1, rhs_psi2, rhs_psi3, rhs_psi5,  vel, Nx);


    //Total concentration initialization
    Vector psi1{Vector::Zero(Nx)};
    Vector psi2{Vector::Zero(Nx)};
    Vector psi3{Vector::Zero(Nx)};
    Vector psi4{Vector::Zero(Nx)};
    Vector psi5{Vector::Zero(Nx)};

    //Reaction term
    Vector rd{Vector::Zero(Nx)};


//Loop Temporale

    for(unsigned int i=1; i<Nt+1; i++)
     {
        concentration.compute_psi(i-1,psi1,psi2,psi3,psi4,psi5); // phi 

        concentration.compute_rd_kp(i-1,rd);//Calcolo i termini di reazione
       
        concentration.one_step_transport_reaction(psi1,psi2,psi3,psi4,psi5,rd,M_rhs,rhs_psi1,rhs_psi2,rhs_psi3,rhs_psi5,i,solver,solver1);
 
        concentration.compute_concentration(i,psi1,psi2,psi3,psi4,psi5); //Calcolo le Concentrazioni effettive
        
    }


    concentration.output_results_fixed_time("Ca");
    concentration.output_results_fixed_time("H_piu");
    concentration.output_results_fixed_time("HCO3_meno");
    concentration.output_results_fixed_time("CO2");
    concentration.output_results_fixed_time("CaSiO3");
    concentration.output_results_fixed_time("SiO2"); 

    //concentration.output_all_reagents(Nx-1);

}
