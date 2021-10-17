#ifndef TYPES_HH
#define TYPES_HH

#include <Eigen/Dense>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> Matrix;

typedef Eigen::MatrixXd Matrix_full;

typedef Eigen::VectorXd Vector;

typedef Eigen::SparseLU<Matrix, Eigen::COLAMDOrdering<int> > Solver;






#endif

