#ifndef SYSTEMS_HH
#define SYSTEMS_HH

//#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "matrix.hpp"
#include "parameters.hpp"


typedef Eigen::SparseMatrix<double> Matrix;
typedef Eigen::Triplet<double> Triplet;



/*!
 *Function for definition of the Darcy System
*/
void set_Darcy_system(const Data_Darcy &data, Eigen::SparseMatrix<double> &Matrix,Eigen::VectorXd &rhs, double h);

/*!
 *Function for the definition of the big sparse block matrix
*/
Matrix block_matrix(const Matrix& A, const Matrix& B);


/*!
 *Auxiliar function for the definiton of the big block matrix
*/
void triplets_with_shift(std::vector<Triplet>& t, const Matrix& A, int shift_row, int shift_col);

#endif

