#include "systems.hpp"
#include <iostream>
#include <vector>
#include <utility>

//BISOGNEREBBE CERCARE COME CONCATENARE MATRICI SPARSE, MA NON È COSÌ SEMPLICE



//void set_Darcy_system(const Data_Darcy &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs, double h)
void set_Darcy_system(const Data_Darcy &data, Eigen::SparseMatrix<double> &M,Eigen::VectorXd &rhs, double h)
{
//All the data that are needed to define the Darcy System are extracted from the data structure
    auto &[L, K, phi, mu, Q_in, Q_out, p_in, p_out, f, Nx, BC_in, BC_out]=data;
   
 

//Computation of the spatial step
    

//Assembling of the Mass Matrix of the first equation with its Boundary Condition
    Matrix_A A(Nx+1,Nx+1);

    A.assemble_matrix(K,h,mu,BC_in,BC_out,p_in,p_out,Q_in,Q_out);


std::cout<<"A definito"<<std::endl;
//Assembling of the B Matrix for the first equation of the system and assembling of the -B^T for the second equation of the system
//For B we have to separate all the steps (we can't use the assemble_matrix() function since for -B^T the modification due to the BC condition have not to be set)
    Matrix_B B(Nx+1,Nx);//è diversa dalle note, la definisco come la B^T delle note
    B.set_data(BC_in,BC_out,f,h);
    B.define_matrix();
    //Eigen::MatrixXd B_T{B.get_matrix().transpose()};
    //Eigen::SparseMatrix<double> B_T{B.get_matrix().transpose()};
    B.set_rhs();
 
   std::cout<<"B definito"<<std::endl;
 
   M=block_matrix(A.get_matrix(),B.get_matrix());

    if(BC_in=="Flow")
            M.coeffRef(0,Nx+1)=0.;
    
    if(BC_out=="Flow")
            M.coeffRef(Nx,Nx+Nx)=0.;



/*//Definition of the all matrix of the Darcy System M=(A,B;-B^T,0)
    //Eigen::MatrixXd A_=A.get_matrix();
    //Eigen::MatrixXd B_=B.get_matrix();
    
//Definition of the first row of M (A,B), we call it C
    Eigen::MatrixXd C(A.get_matrix().rows(),A.get_matrix().cols()+B.get_matrix().cols());
    //Eigen::SparseMatrix<double> C(A.get_matrix().rows(),A.get_matrix().cols()+B.get_matrix().cols());
    C<<A.get_matrix(),B.get_matrix();

std::cout<<"C definito"<<std::endl;

//Definition of the second row of M (-B^T,0), we call it E
    //Eigen::MatrixXd D= Eigen::MatrixXd::Zero(B_T.rows(),B.get_matrix().cols());
   Eigen::MatrixXd E(B_T.rows(),B_T.cols()+B_T.rows());
   //Eigen::SparseMatrix<double> E(B_T.rows(),B_T.cols()+B_T.rows());
   E<<-B_T,Eigen::MatrixXd::Zero(B_T.rows(),B.get_matrix().cols());

std::cout<<"E definito"<<std::endl;

//Assembling of the matrix M
    Eigen::MatrixXd M(C.rows()+E.rows(),C.cols());
   
    M<<std::move(C),std::move(E);

 
    //Matrix<<C,E;
std::cout<<"M definito"<<std::endl;*/

//Definition of the rhs of the Darcy System. The rhs of the first equation and the one of the second are glued to build the all rhs.
    //Eigen::VectorXd v1=A.get_rhs();
    //Eigen::VectorXd v2=B.get_rhs();
    //Eigen::VectorXd v3{A.get_rhs().size()+B.get_rhs().size()};
    rhs<<A.get_rhs(),B.get_rhs();


}




void triplets_with_shift(std::vector<Triplet>& t, const Matrix& A, int shift_row, int shift_col)
{
    for (int k(0); k < A.outerSize(); ++k)
        for (Matrix::InnerIterator it(A,k); it; ++it)
            t.push_back(Triplet(it.row() + shift_row, it.col() + shift_col, it.value()));
}

Matrix block_matrix(const Matrix& A, const Matrix& B)
{
    // M = [[A, B]; [B^T, 0]]
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




