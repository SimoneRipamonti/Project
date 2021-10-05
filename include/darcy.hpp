#ifndef SYSTEMS_HH
#define SYSTEMS_HH

//#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "matrix.hpp"
#include "parameters.hpp"



typedef Eigen::Triplet<double> Triplet;


/** \addtogroup Darcy_Functions
 *\brief Functions that are related to the definition and resolution of the Darcy system
 *  @{
 */
/*!
 *Function for definition of the Darcy System
*\param data is a reference to the data that we need to define the Darcy problem
*\param Matrix is a reference to a the sparse matrix where we want to store the entire Darcy system matrix ([[A,B];[-B^T,0]])
*\param rhs is a reference to the rhs of the Darcy system
*\param h is the spatial step
*/
void set_Darcy_system(const Data_Darcy &data, Matrix &Matrix, Vector &rhs, double h);

/*!
 *Function for the definition of the big sparse block matrix
*\param A is the A matrix of the Darcy system
*\param B is the B matrix of the Darcy system
*\param Matrix is the output darcy system matrix
*/
Matrix block_matrix(const Matrix& A, const Matrix& B);


/*!
 *Auxiliar function for the definiton of the big block matrix
*/
void triplets_with_shift(std::vector<Triplet>& t, const Matrix& A, int shift_row, int shift_col);

/*!
 *Function that creates and solves the Darcy system in order to compute the Darcy velocity (it is used in case6_all) 
*\param data is a reference to the data that we need to define the Darcy problem
*\param vel is a reference to the vector where we want to store the Darcy velocity
 */
void Darcy_velocity(const Data_Darcy &data, Vector &vel);

/** @}*/









#endif

