#ifndef CONCENTRATIONS_HH
#define CONCENTRATIONS_HH

#include <Eigen/Dense>
#include "parameters.hpp"

#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <string>


enum Method{EsplicitEuler,PredictorCorrector,Heun};//Numerical scheme for solving the reaction part

class Concentration
{
public:
   Concentration(const std::string &filename);

   unsigned int get_Nx();
  
   unsigned int get_Nt();

   void set_initial_cond(); //set initial condition for the reagents

   void assemble_transport(Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& vel); //set the matrix for the transport part of the equation
 
   void  assemble_transport_CO2(Eigen::SparseMatrix<double>& M_CO2, Eigen::VectorXd& rhs_CO2, const Eigen::VectorXd& vel, double C_in, double C_out, const std::string& bc_cond); //set the matrix for the CO2 transport part of the equation

   void compute_phi(unsigned int step, Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5); //compute phi variables from reagents

   void compute_rd(unsigned int step, Eigen::VectorXd& rd);//compute the reaction term of the equation

   void compute_rd_kp(unsigned int step, Eigen::VectorXd& rd);//compute the reaction term of the equation

   void one_step_transport_reaction(Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5, Eigen::VectorXd& rd, const  Eigen::SparseMatrix<double> rhs, const Eigen::VectorXd& rhs_CO2, unsigned int step, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2); //compute one temporal step for the transport-reaction equation

   void Euler_Esplicit(Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5, const Eigen::VectorXd& rd, const Eigen::SparseMatrix<double>&  rhs, const Eigen::VectorXd& rhs_CO2, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2); //solve the reaction part with an Esplicit Euler

   void transport_and_reaction(Eigen::VectorXd& psi, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rd,Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver);
   //solve the equation for one psi

   void transport_and_reaction_CO2(Eigen::VectorXd& psi, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rhs_CO2, const Eigen::VectorXd& rd, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2);//solve the equation for the total concentration that include CO2 (psi2)

   //Risoluzione sistema non lineare
   
   void compute_concentration(unsigned int step, const Eigen::VectorXd& psi1, const Eigen::VectorXd& psi2, const Eigen::VectorXd& psi3, const Eigen::VectorXd& psi4, const Eigen::VectorXd& psi5); //Solve the non linear system with the Newthon method, in order to compute the real concentrations

   void compute_rhs(Eigen::VectorXd& rhs, const Eigen::VectorXd& old_it, double phi1, double phi2, double phi3, double phi4, double phi5); //Compute the rhs of the Newton Scheme

   void compute_Jacob(Eigen::MatrixXd& J,const Eigen::VectorXd& old_it); //Compute the unknowns Jacobian Element of the non linear system

   //These functions save the results in .csv files

   void output_results_fixed_time(const std::string&); //row=space, column=time

   void output_results_fixed_space(const std::string&);//row=time, column=space

   void output_all_reagents(unsigned int pos);//Evolution in time of the species in a given position


private:
   Eigen::MatrixXd Ca;
   Eigen::MatrixXd H_piu;
   Eigen::MatrixXd HCO3_meno;
   Eigen::MatrixXd CO2;
   Eigen::MatrixXd CaSiO3;
   Eigen::MatrixXd SiO2;

   Data_Transport data_transp;
   Data_6Reagents data_reagents;
   Data_Reaction data_reaction;

   Method method;
};







#endif
