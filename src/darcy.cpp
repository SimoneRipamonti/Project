#include "darcy.hpp"
#include <iostream>
#include <vector>
#include <utility>

void set_Darcy_system(const Data_Darcy &data, Matrix &M, Vector &rhs, double h)
{
//All the data that are needed to define the Darcy System are extracted from the data structure
    auto &[L, K, phi, mu, Q_in, Q_out, p_in, p_out, f, Nx, BC_in, BC_out]=data;


//Assembling of the Mass Matrix of the first equation with its Boundary Condition
    Matrix_A A(Nx+1,Nx+1);
    A.assemble_matrix(K,h,mu,BC_in,BC_out,p_in,p_out,Q_in,Q_out);

//Assembling of the B matrix of the Darcy system
    Matrix_B B(Nx+1,Nx);
    B.set_data(BC_in,BC_out,f,h);
    B.define_matrix();
    B.set_rhs();

//Definition of the Darcy system big sparse matrix
    M=block_matrix(A.get_matrix(),B.get_matrix());

//Setting of the BC for the big matrix
    if(BC_in=="Flow")
        M.coeffRef(0,Nx+1)=0.;

    if(BC_out=="Flow")
        M.coeffRef(Nx,Nx+Nx)=0.;

//Definition of the rhs of the Darcy system
    rhs<<A.get_rhs(),B.get_rhs();


}


Matrix block_matrix(const Matrix& A, const Matrix& B)
{
    std::vector<Triplet> t;
    t.reserve(A.nonZeros() + B.nonZeros());

    triplets_with_shift(t, A, 0, 0);
    triplets_with_shift(t, B, 0, A.cols());
    triplets_with_shift(t, -B.transpose(), A.rows(), 0);

    Matrix M(A.rows() + B.cols(), A.cols() + B.cols());
    M.setFromTriplets(std::begin(t), std::end(t));
    M.makeCompressed();

    return M;
}


void triplets_with_shift(std::vector<Triplet>& t, const Matrix& A, int shift_row, int shift_col)
{
    for (int k(0); k < A.outerSize(); ++k)
        for (Matrix::InnerIterator it(A,k); it; ++it)
            t.push_back(Triplet(it.row() + shift_row, it.col() + shift_col, it.value()));
}




void Darcy_velocity(const Data_Darcy &data, Vector &vel)
{
  
    double h =static_cast<double>(data.L)/data.Nx; //space step

    Eigen::VectorXd sol(data.Nx+data.Nx+1); //solution vectors are resized

    Eigen::SparseMatrix<double> M_d(data.Nx+data.Nx+1,data.Nx+data.Nx+1);//Initialization of the big matrix for the Darcy system
    Eigen::VectorXd rhs_d(data.Nx+data.Nx+1);//Initialization of the rhs of Darcy
    set_Darcy_system(data,M_d,rhs_d,h);//Definition of the Darcy system Mx=rhs


    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver_d; //Initialization of the Solver for the sparse system

    solver_d.analyzePattern(M_d); // Compute the ordering permutation vector from the structural pattern of A

    solver_d.factorize(M_d); // Compute the numerical factorization

    sol= solver_d.solve(rhs_d);//The Darcy system is solved and the solution is stored in the sol vector

    vel=sol.head(data.Nx+1);


}




