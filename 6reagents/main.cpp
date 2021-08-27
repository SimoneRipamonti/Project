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



int main()
{

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

    Eigen::VectorXd vel{0.0*Eigen::VectorXd::Ones(Nx+1)};


    Eigen::VectorXd psi1{Eigen::VectorXd::Zero(Nx)};
    Eigen::VectorXd psi2{Eigen::VectorXd::Zero(Nx)};
    Eigen::VectorXd psi3{Eigen::VectorXd::Zero(Nx)};
    Eigen::VectorXd psi4{Eigen::VectorXd::Zero(Nx)};
    Eigen::VectorXd psi5{Eigen::VectorXd::Zero(Nx)};



//Set the initial condition

    concentration.set_initial_cond();


//Setting for transport part of the equation
    Eigen::SparseMatrix<double> M(Nx,Nx);//Mass matrix for the transport part
    Eigen::SparseMatrix<double> rhs(Nx,Nx);//rhs for the transport part
    concentration.assemble_transport(M,rhs,vel);//function that assembles the transport part


//CO2 has to be treat separately since it is added constantly in the input section
    auto [C_in,C_out,bc_cond]=data_CO2;
    Eigen::SparseMatrix<double> M_CO2(Nx,Nx);
    Eigen::VectorXd rhs_CO2(Nx);
    concentration.assemble_transport_CO2(M_CO2,rhs_CO2,vel,C_in,C_out,bc_cond);


//Here solver for the transport system is initialized
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    solver.analyzePattern(M);
    solver.factorize(M);

//We make the same also for the CO2 transport part
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver2;
    solver2.analyzePattern(M_CO2);
    solver2.factorize(M_CO2);

//Reaction term
    Eigen::VectorXd rd{Eigen::VectorXd::Zero(Nx)};


//Loop Temporale

    for(unsigned int i=1; i<Nt; i++)
    {
        concentration.compute_phi(i-1,psi1,psi2,psi3,psi4,psi5); // phi

        concentration.compute_rd(i-1,rd);//Calcolo i termini di reazione

        //concentration.one_step_transport_reaction(phi1,phi2,phi3,phi4,phi5,rd,M,rhs,method,i,solver); //Calcolo un passo della Reazione

        concentration.one_step_transport_reaction(psi1,psi2,psi3,psi4,psi5,rd,rhs,rhs_CO2,i,solver,solver2);

        concentration.compute_concentration(i,psi1,psi2,psi3,psi4,psi5); //Calcolo le Concentrazioni effettive
    }


    concentration.output_results_fixed_time("Ca");

    concentration.output_results_fixed_space("Ca");

    concentration.output_results_fixed_space("CaSiO3");

    concentration.output_all_reagents(Nx-1);

}
