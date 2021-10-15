#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"



void set_solver(Matrix& M, Solver& solver);
/*!
 *Definition of the esplicit tranport system 
 */
void Transport_esp(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond, Matrix &M, Matrix &rhs1, Vector &rhs2);

/*!
 *Definition of the Implicit transport system 
 */
void Transport_imp(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond, Matrix &M, Matrix &rhs1, Vector &rhs2);

#endif

