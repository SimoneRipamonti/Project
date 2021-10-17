#ifndef FUNCTIONS_HH
#define FUNCTIONS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

/** \addtogroup Functions
 *\brief Auxiliar Functions used in the code
 *  @{
*/
/*!
 *Function that sets the solver for the solution of the linear system denoted by the matrix M
 *\param M is the sparse matrix representing the linear system that has to be solved
 *\param solver solver for the linear systems that uses the LU factorization 
*/
void set_solver(Matrix& M, Solver& solver);

/*!
 *Function that sets the explicit transport system 
 */
void Transport_exp(Matrix_full &Ca, Vector &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond, Matrix &M, Matrix &rhs1, Vector &rhs2);

/*!
 *Function that sets the implicit transport system
 */
void Transport_imp(Matrix_full &Ca, Vector &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond, Matrix &M, Matrix &rhs1, Vector &rhs2);

/** @}*/

#endif

