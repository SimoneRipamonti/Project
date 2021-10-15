#include "functions.hpp"
#include "transport_decay.hpp"
#include <iostream>
#include <vector>


//Esplicit transport upwind
void Transport_system_esplicit(Eigen::MatrixXd &Ca, Vector &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond) //As input there is the matrix solution where we store our solution at each istant: each row represent a spatial position, each column represent a time istant. In the vector vel there is the velocity evaluated at each node cell.
{
    unsigned int Nx{data_transport.Nx};
    Matrix M(Nx,Nx), rhs1(Nx,Nx);
    Vector rhs2(Nx);
    //Definition of the linear system explicit matrix
    Transport_esp(Ca,vel,data_transport,initial_cond,M,rhs1,rhs2);  
    
    Solver solver;
    set_solver(M,solver);
    
    //Initialization of the rhs of the Transport System
    Vector rhs(Nx);

    //Temporal loop for solving at each istant the transport problem
    for(unsigned int i=1; i<data_transport.Nt+1; i++)
    {
        rhs=rhs1*Ca.col(i-1)+rhs2;
        Ca.col(i)=solver.solve(rhs);
    }
}






//Implicit transport upwind and esplicit reaction
void Transport_system_implicit(Eigen::MatrixXd &Ca, Vector &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond) //As input there is the matrix solution where we store our solution at each istant,
//each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{
    unsigned int Nx{data_transport.Nx};
    Matrix M(Nx,Nx), rhs1(Nx,Nx);
    Vector rhs2(Nx);
    //Definition of the linear system explicit matrix
    Transport_imp(Ca,vel,data_transport,initial_cond,M,rhs1,rhs2);  
    
    Solver solver;
    set_solver(M,solver);

    //Initialization of the rhs of the Transport System
    Vector rhs(Nx);

    //Temporal loop for solving at each istant the transport problem
    for(unsigned int i=1; i<data_transport.Nt+1; i++)
    {
        rhs=rhs1*Ca.col(i-1)+rhs2;
        Ca.col(i)=solver.solve(rhs);
    }
}





