#ifndef SYSTEMS_HH
#define SYSTEMS_HH

//#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "matrix.hpp"
#include "parameters.hpp"


typedef Eigen::SparseMatrix<double> Matrix;
typedef Eigen::Triplet<double> Triplet;


//Definition of the Darcy System
//void set_Darcy_system(const Data_Darcy &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs, double h);
void set_Darcy_system(const Data_Darcy &data, Eigen::SparseMatrix<double> &Matrix,Eigen::VectorXd &rhs, double h);

void triplets_with_shift(std::vector<Triplet>& t, const Matrix& A, int shift_row, int shift_col);

Matrix block_matrix(const Matrix& A, const Matrix& B);
/*//Definition of the Transport system solved with an explicit scheme
void Transport_system_esplicit(Eigen::MatrixXd &solution,Eigen::VectorXd &velocity,Data_Transport &data);

//Definition of the Tranport system solved with an implicit scheme
void Transport_system_implicit(Eigen::MatrixXd &Ca,Eigen::MatrixXd &CaSiO3,Eigen::VectorXd &velocity,Data_Transport &data_transport, Data_Reaction &data_reaction);
*/
#endif

