#include <Eigen/Dense>

void set_initial_cond(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, Eigen::MatrixXd& SiO2, muparser_fun Ca0, muparser_fun H_piu0,
                      muparser_fun HCO3_meno0, muparser_fun CO20, muparser_fun CaSiO30, double h);

void assemble_transport(Eigen::MatrixXd& M, Eigen::VectorXd& rhs, Eigen::VectorXd vel, muparser_fun phi, double h, unsigned int Nx);

void compute_phi(Eigen::VectorXd& phi1,Eigen::VectorXd& phi2,Eigen::VectorXd& phi3,Eigen::VectorXd& phi4,Eigen::VectorXd& phi5,Eigen::VectorXd Ca,Eigen::VectorXd H_piu,Eigen::VectorXd HCO3_meno, Eigen::VectorXd CO2,Eigen::VectorXd CaSiO3,Eigen::VectorXd SiO2);

void compute_rd(&rd, const Eigen::VectorXd Ca, const Eigen::VectorXd H_piu, const Eigen::VectorXd SiO2, double const_r, double K_eq, double n);

void one_step_transport_reaction(Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, Eigen::VectorXd& phi5, const Eigen::VectorXd rd, const auto M_lu, Eigen::MatrixXd rhs);

void transport_and_reaction(Eigen::VectorXd& phi, const auto M_lu, const Eigen::MatrixXd &rhs, const Eigen::VectorXd rd);

//Risoluzione sistema non lineare
void compute_concentration(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, Eigen::MatrixXd& SiO2, const Eigen::VectorXd phi1, const Eigen::VectorXd phi2, const Eigen::VectorXd phi3, const Eigen::VectorXd phi4, const Eigen::VectorXd phi5,unsigned int step,double K_sol);

void compute_rhs(Eigen::VectorXd& rhs,const Eigen::VectorXd old_it, double phi1, double phi2, double phi3, double phi4, double phi5, double K_eq);

void compute_Jacob_last_row(Eigen::MatrixXd& J,const Eigen::VectorXd old_it);

