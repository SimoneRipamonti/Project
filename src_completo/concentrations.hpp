#ifndef CONCENTRATIONS_HH
#define CONCENTRATIONS_HH

#include <Eigen/Dense>
#include "parameters.hpp"

#include <string>
class Concentration
{
public:
   Concentration(const std::string &filename);

   unsigned int get_Nx();
  
   unsigned int get_Nt();

   /*double get_L();
 
   double get_T();*/

   void set_initial_cond();

   void assemble_transport(Eigen::MatrixXd& M, Eigen::MatrixXd& rhs, const Eigen::VectorXd& vel);

   void compute_phi(unsigned int step, Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, Eigen::VectorXd& phi5);

   void compute_rd(unsigned int step, Eigen::VectorXd& rd);

   












void one_step_transport_reaction(Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, Eigen::VectorXd& phi5, Eigen::VectorXd& rd, const Eigen::MatrixXd& M, const Eigen::MatrixXd& rhs, int method, unsigned int step);
   
void Euler_Esplicit(Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, Eigen::VectorXd& phi5, const Eigen::VectorXd& rd, const Eigen::MatrixXd& M, const Eigen::MatrixXd& rhs);

void transport_and_reaction(Eigen::VectorXd& phi, const Eigen::MatrixXd& M, const Eigen::MatrixXd& rhs, const Eigen::VectorXd& rd);

   //Risoluzione sistema non lineare
   

























   void compute_concentration(unsigned int step, const Eigen::VectorXd& phi1, const Eigen::VectorXd& phi2, const Eigen::VectorXd& phi3, const Eigen::VectorXd& phi4, const Eigen::VectorXd& phi5);

   void compute_rhs(Eigen::VectorXd& rhs, const Eigen::VectorXd& old_it, double phi1, double phi2, double phi3, double phi4, double phi5);

   void compute_Jacob(Eigen::MatrixXd& J,const Eigen::VectorXd& old_it);

   void output_results_fixed_time(const std::string&);

   void output_results_fixed_space(const std::string&);

   void output_all_reagents(unsigned int pos);


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
};







#endif
