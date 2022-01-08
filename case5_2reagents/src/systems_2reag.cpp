#include "systems_2reag.hpp"
#include <iostream>
#include <vector>

//Explicit transport upwind and explicit reaction
void Transport_system_esplicit_2_reagents(Matrix_full &Ca, Matrix_full &CaSiO3, Vector &vel,Data_Transport &data_transport, Data_Reaction &data_reaction, Data_2Reagents &data_2reagents) //As input there is the matrix solution where we store our solution at each istant,
//each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,method]=data_transport;
    auto &[Ca_0,CaSiO3_0,K_sol,ph]=data_2reagents;
    auto &[Area,rate_const,E,R,Temperature]=data_reaction;

    //Computation of the spatial and temporal step from the data
    const double h =static_cast<double>(L)/Nx;
    const double dt=static_cast<double>(T)/Nt;


    //The Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System)

    for (unsigned int i=0; i<Nx; ++i)
        Ca(i,0)=Ca_0(h/2+i*h);

    //Ca(0,0)=1.0;
    for (unsigned int i=0; i<Nx; ++i)
        CaSiO3(i,0)=CaSiO3_0(h/2+i*h);

    //Definition of the upwind transport matrices
    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    //Definition of the mass transport matrix
    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);

    //Definition of the Reaction Matrix
    Matrix_R React(Nx,Nx);
    React.assemble_matrix(Area,rate_const,Temperature,R,E,ph,K_sol,h,phi);

    Matrix M{1/dt*C.get_matrix()};//matrix of the transport linear system

    Solver  solver; //solver used to solve the sparse linear system
 
    set_solver(M,solver);

    Vector rhs(Nx);//rhs of the transport linear system

    //Temporal loop
    for(unsigned int i=1; i<Nt+1; i++)
    {
        React.update(Ca.col(i-1));
        rhs=(1/dt*C.get_matrix()-F_p.get_matrix()+F_m.get_matrix())*Ca.col(i-1)+F_p.get_rhs()-F_m.get_rhs()+React.get_rhs().cwiseProduct(CaSiO3.col(i-1));;
        Ca.col(i)=solver.solve(rhs);
        CaSiO3.col(i)=CaSiO3.col(i-1)-dt/(h*phi(h/2+i*h))*React.get_rhs().cwiseProduct(CaSiO3.col(i-1));
    }
}









//Implicit transport upwind and explicit reaction
void Transport_system_implicit_2_reagents(Matrix_full &Ca, Matrix_full &CaSiO3, Vector &vel,Data_Transport &data_transport,Data_Reaction &data_reaction, Data_2Reagents &data_2reagents) //As input there is the matrix solution where we store our solution at each istant,
//each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond,method]=data_transport;
    auto &[Ca_0,CaSiO3_0,K_sol,ph]=data_2reagents;
    auto &[Area,rate_const,E,R,Temperature]=data_reaction;

    //Computation of the spatial and temporal step from the data
    const double h =static_cast<double>(L)/Nx;
    const double dt=static_cast<double>(T)/Nt;


    //The Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System)

    for (unsigned int i=0; i<Nx; ++i)
        Ca(i,0)=Ca_0(h/2+i*h);

    for (unsigned int i=0; i<Nx; ++i)
        CaSiO3(i,0)=CaSiO3_0(h/2+i*h);

    //Definition of the upwind transport matrices
    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    //Definition of the mass transport matrix
    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);

    //Definition of the reaction matrix
    Matrix_R React(Nx,Nx);
    React.assemble_matrix(Area,rate_const,Temperature,R,E,ph,K_sol,h,phi);


    Matrix M{1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix()};//matrix of the transport linear system

    Solver   solver; //solver used to solve the sparse linear system

    set_solver(M, solver);
    
    Vector rhs(Nx);//rhs of the transport linear system

    //Temporal loop
    for(unsigned int i=1; i<Nt+1; i++)
    {
        React.update(Ca.col(i-1));
        rhs=(1/dt*C.get_matrix())*Ca.col(i-1)+F_p.get_rhs()-F_m.get_rhs()+React.get_rhs().cwiseProduct(CaSiO3.col(i-1));
        Ca.col(i)=solver.solve(rhs);
        CaSiO3.col(i)=CaSiO3.col(i-1)-dt/(h*phi(h/2+i*h))*React.get_rhs().cwiseProduct(CaSiO3.col(i-1));
    }
}





