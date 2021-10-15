#include "functions.hpp"
#include "transport_decay.hpp"
#include <iostream>
#include <vector>


//Explicit transport upwind and explicit reaction
void Transport_system_explicit(Eigen::MatrixXd &Ca, Vector &vel, Data_Transport &data_transport, Data_linear_decay &data_linear_decay)
{
    //Initialization of the matrices and vectors for the transport and decay problem
    unsigned int Nx{data_transport.Nx};
    Matrix M(Nx,Nx), rhs1(Nx,Nx);
    Vector rhs2(Nx);
    
    //Definition of the matrices and vector needed to solve the transport and decay problem
    Transport_exp(Ca,vel,data_transport,data_linear_decay,M,rhs1,rhs2);  
    
    //Set of the linear system solver
    Solver solver;
    set_solver(M,solver);
    
    //Initialization of the rhs of the transport and decay system
    Vector rhs(Nx);

    //Temporal loop for solving at each istant the transport problem
    for(unsigned int i=1; i<data_transport.Nt+1; i++)
    {
        rhs=rhs1*Ca.col(i-1)+rhs2;
        Ca.col(i)=solver.solve(rhs);
    }
}






//Implicit transport upwind and esplicit reaction
void Transport_system_implicit(Eigen::MatrixXd &Ca, Vector &vel, Data_Transport &data_transport, Data_linear_decay &data_linear_decay) 
{
    //Initialization of the matrices and vectors for the transport and decay problem
    unsigned int Nx{data_transport.Nx};
    Matrix M(Nx,Nx), rhs1(Nx,Nx);
    Vector rhs2(Nx);
    
    //Definition of the matrices and vector needed to solve the transport and decay problem
    Transport_imp(Ca,vel,data_transport,data_linear_decay,M,rhs1,rhs2);  
    
    //Set of the linear system solver
    Solver solver;
    set_solver(M,solver);

    //Initialization of the rhs of the transport and decay system
    Vector rhs(Nx);

    //Temporal loop for solving at each istant the transport problem
    for(unsigned int i=1; i<data_transport.Nt+1; i++)
    {
        rhs=rhs1*Ca.col(i-1)+rhs2;
        Ca.col(i)=solver.solve(rhs);
    }
}





