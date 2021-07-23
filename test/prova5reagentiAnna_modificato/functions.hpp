#include <Eigen/Dense>
#include "muparser_fun.hpp"
#include "matrix.hpp"

#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

void change_element (Eigen::MatrixXd& a);

void set_initial_cond(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, muparser_fun Ca0, muparser_fun H_piu0,
                      muparser_fun HCO3_meno0, muparser_fun CO20, muparser_fun CaSiO30, double h);

void assemble_transport(Eigen::SparseMatrix<double>& M, Eigen::MatrixXd& rhs, Eigen::VectorXd vel, muparser_fun phi, double h, unsigned int Nx, double dt);

void compute_phi(Eigen::VectorXd& phi1,Eigen::VectorXd& phi2,Eigen::VectorXd& phi3,Eigen::VectorXd& phi4, Eigen::VectorXd Ca,Eigen::VectorXd H_piu,Eigen::VectorXd HCO3_meno, Eigen::VectorXd CO2,Eigen::VectorXd CaSiO3);

void compute_rd(Eigen::VectorXd& rd, const Eigen::VectorXd Ca, const Eigen::VectorXd H_piu, double const_r, double K_eq, double n);

void one_step_transport_reaction(Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, const Eigen::VectorXd rd, const Eigen::MatrixXd rhs, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver);

void transport_and_reaction(Eigen::VectorXd& phi, const Eigen::MatrixXd rhs, const Eigen::VectorXd rd, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver);

//Risoluzione sistema non lineare
void compute_concentration(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, const Eigen::VectorXd phi1, const Eigen::VectorXd phi2, const Eigen::VectorXd phi3, const Eigen::VectorXd phi4, unsigned int step,double K_eq);

void compute_rhs(Eigen::VectorXd& rhs,const Eigen::VectorXd old_it, double phi1, double phi2, double phi3, double phi4, double K_eq);

void compute_Jacob(Eigen::MatrixXd& J,const Eigen::VectorXd old_it);

