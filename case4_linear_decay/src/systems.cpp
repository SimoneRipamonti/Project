#include "systems.hpp"
#include <iostream>
#include <vector>

#include <Eigen/LU>

#include <Eigen/Sparse>
#include <Eigen/SparseLU>

//Esplicit transport upwind

void Transport_system_esplicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel,Data_Transport &data_transport, Data_linear_decay &data_linear_decay) //As input there is the matrix solution where we store our solution at each istant, each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,method]=data_transport;
    auto &[Ca_0,lambda]=data_linear_decay;

    //Computation of the spatial and temporal step from the data
    double h =static_cast<double>(L)/Nx;
    double dt=static_cast<double>(T)/Nt;

    Eigen::MatrixXd Reaction{h*Eigen::MatrixXd::Identity(Nx,Nx)}; //Reaction matrix for the linear decay (It is just the Identity matrix multiplied by h and the porosity)

    //The Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System) 
    for (unsigned int i=0; i<Nx; ++i)
        {Ca(i,0)=Ca_0(h/2+i*h);
         Reaction(i,i)*=phi(h/2+i*h);}

    //Definition of the Upwind matrices
    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    //Definition of the Mass Matrix for the transport
    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);

    //Definition of the linear system explicit matrix
    const Eigen::SparseMatrix<double>  M{1/dt*C.get_matrix()}; //Definition of the sparse upwind matrix for the esplicit scheme

    //Initialization of the rhs of the Transport System
    Eigen::VectorXd rhs(Nx);//Initialization of the rhs of the system

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver; //Solver fot the sparse system

    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(M);
    // Compute the numerical factorization
    solver.factorize(M);

     //Temporal loop for solving at each istant the linear decay\transport problem
    for(unsigned int i=1; i<Nt+1; i++)
    {
        rhs=(1/dt*C.get_matrix()-F_p.get_matrix()+F_m.get_matrix()-lambda*Reaction)*Ca.col(i-1)+C.get_rhs();
        Ca.col(i)=solver.solve(rhs);

    }

}

//Implicit transport upwind and esplicit reaction
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel,Data_Transport &data_transport, Data_linear_decay &data_linear_decay) //As input there is the matrix solution where we store our solution at each istant, each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,method]=data_transport;
    auto &[Ca_0,lambda]=data_linear_decay;

    //Computation of the spatial and temporal step from the data
    double h =static_cast<double>(L)/Nx;
    double dt=static_cast<double>(T)/Nt;


    Eigen::MatrixXd Reaction{h*Eigen::MatrixXd::Identity(Nx,Nx)}; //Reaction matrix for the linear decay (It is just the Identity matrix multiplied by h and the porosity)

   //Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System) 
   for (unsigned int i=0; i<Nx; ++i)
        {Ca(i,0)=Ca_0(h/2+i*h);
         Reaction(i,i)*=phi(h/2+i*h);}

    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);

    const Eigen::SparseMatrix<double>  M{1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix()}; //Definition of the sparse upwind matrix for the implicit scheme

    //Initialization of the rhs of the system
    Eigen::VectorXd rhs(Nx);

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver; //Solver fot the sparse system

    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(M);
    // Compute the numerical factorization
    solver.factorize(M);

    //Temporal loop for solving at each istant the linear decay\transport problem
    for(unsigned int i=1; i<Nt+1; i++)
    {
        rhs=(1/dt*C.get_matrix()-lambda*Reaction)*Ca.col(i-1)+C.get_rhs();
        Ca.col(i)=solver.solve(rhs);

    }

}



