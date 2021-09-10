#include "systems.hpp"
#include <iostream>
#include <vector>


//Esplicit transport upwind
void Transport_system_esplicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel,Data_Transport &data_transport, Data_linear_decay &initial_cond) //As input there is the matrix solution where we store our solution at each istant: each row represent a spatial position, each column represent a time istant. In the vector vel there is the velocity evaluated at each node cell.
{
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,method]=data_transport;

    auto &[C_0,lambda]=initial_cond;    
 
    //Computation of the spatial and temporal step from the data
    const double h =static_cast<double>(L)/Nx;
    const double dt=static_cast<double>(T)/Nt;


    //The Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System)

    for (unsigned int i=0; i<Nx; ++i)
        Ca(i,0)=C_0(h/2+i*h);
    
    //We need to have compatibility between the bc and temproal ones for the first concentration degree of freedom
    if(bc_cond=="In")
       Ca(0,0)=C_in;
    else
       Ca(Nx-1,Nt-1)=C_out;


    //Definition of the Upwind matrices
    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    //Definition of the Mass Matrix for the transport
    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);

    //Definition of the linear system explicit matrix
    const Eigen::SparseMatrix<double> M{1/dt*C.get_matrix()};

    //Solver for the sparse Transport system
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(M);
    // Compute the numerical factorization
    solver.factorize(M);

    //Initialization of the rhs of the Transport System
    Eigen::VectorXd rhs(Nx);

    //Temporal loop for solving at each istant the transport problem
    for(unsigned int i=1; i<Nt+1; i++)
    {
        rhs=(1/dt*C.get_matrix()-F_p.get_matrix()+F_m.get_matrix())*Ca.col(i-1)+F_p.get_rhs()-F_m.get_rhs();
        Ca.col(i)=solver.solve(rhs);
    }
}






//Implicit transport upwind and esplicit reaction
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond) //As input there is the matrix solution where we store our solution at each istant,
//each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,meth]=data_transport;
    //muparser_fun C_0{initial_cond.C_0};
     auto &[C_0,lambda]=initial_cond; 

    //Computation of the spatial and temporal step from the data
    const double h =static_cast<double>(L)/Nx;
    const double dt=static_cast<double>(T)/Nt;


    //The Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System)

    for (unsigned int i=0; i<Nx; ++i)
        Ca(i,0)=C_0(h/2+i*h);
   
    if(bc_cond=="In")
       Ca(0,0)=C_in;
    else
       Ca(Nx-1,Nt-1)=C_out;


    //Definition of the Upwind matrices
    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    //Definition of the Mass Matrix for the transport
    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);

    //Definition of the linear system implicit matrix
    const Eigen::SparseMatrix<double> M{1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix()};

    //Solver for the sparse Transport system
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(M);
    // Compute the numerical factorization
    solver.factorize(M);

    //Initialization of the rhs of the Transport System
    Eigen::VectorXd rhs(Nx);

    //Temporal loop for solving at each istant the transport problem
    for(unsigned int i=1; i<Nt+1; i++)
    {
        rhs=(1/dt*C.get_matrix())*Ca.col(i-1)+F_p.get_rhs()-F_m.get_rhs();
        Ca.col(i)=solver.solve(rhs);
    }
}




