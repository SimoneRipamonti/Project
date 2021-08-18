#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

//Function that saves the solution vector of the Darcy system on a csv file
//void Darcy_output_results(Eigen::VectorXd &value,unsigned int Nx,double L);

//Function that saves the solution matrix (each row represent a spatial position and each column a temporal one) on a csv file 
void output_results_fixed_time(const std::string& title, const Eigen::MatrixXd &value1, unsigned int Nx, double L,unsigned int Nt);

void output_results_fixed_space(const std::string& title, const Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);

/*void output_all_reagents(const Eigen::MatrixXd &Ca,const Eigen::MatrixXd &H_piu,const Eigen::MatrixXd &CaSiO3,const Eigen::MatrixXd &CO2,const Eigen::MatrixXd &SiO2,const Eigen::MatrixXd &HCO3_meno, unsigned int position, double T, unsigned int Nt); 

void pressure_exact_result(Eigen::VectorXd &value1, unsigned int Nx, double L);

void velocity_exact_result(Eigen::VectorXd &value1, unsigned int Nx, double L);

void output_error(Eigen::VectorXd &value1, Eigen::VectorXi &Nx);

void Transport_output_results_fixed_time(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2, unsigned int Nx, double L,unsigned int Nt);

void Transport_output_results_fixed_space(Eigen::MatrixXd &value1, Eigen::MatrixXd &value2,unsigned int Nx, double T,unsigned int Nt);
*/


//void Transport_output_results_fixed_space(Eigen::MatrixXd &value1, unsigned int Nx, double T,unsigned int Nt);



#endif

