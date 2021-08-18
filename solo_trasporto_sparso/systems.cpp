#include "systems.hpp"
#include <iostream>
#include <vector>


//Esplicit transport upwind
void Transport_system_esplicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel,Data_Transport &data_transport,Data_2Reagents &data_2reagents) //As input there is the matrix solution where we store our solution at each istant,
//each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{ 
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond]=data_transport;
    auto &[Ca_0,CaSiO3_0,K_sol,ph]=data_2reagents;

    //Computation of the spatial and temporal step from the data
    double h =static_cast<double>(L)/Nx;//constexpr?
    double dt=static_cast<double>(T)/Nt;//constexpr?

    //Porosità posta uguale a 1 
    //muparser_fun phi;
    //std::string s="1.0+0.0*x";
    //phi.set_value(s);

    //The Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System)
   
    for (unsigned int i=0; i<Nx; ++i)
        Ca(i,0)=Ca_0(h/2+i*h);



    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);
    

   //std::cout<<phi(1)<<std::endl;


    
    //Eigen::MatrixXd M(Nx,Nx);//Matrix of the linear system that has to be solved
    //Eigen::VectorXd rhs(Nx);//rhs of the lineay system that has to be solved

    //Here there is a different definition of the linear matrix of the system (and also of the rhs in the time loop)
    const Eigen::SparseMatrix<double> M(1/dt*C.get_matrix());
    //const Eigen::MatrixXd M(1/dt*C.get_matrix());
  
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    // fill A and b;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(M); 
    // Compute the numerical factorization 
    solver.factorize(M); 

    Eigen::VectorXd rhs(Nx);

    //std::cout<<F_p.get_matrix()<<std::endl;
    //std::cout<<F_m.get_matrix()<<std::endl;
    //std::cout<<C.get_matrix()<<std::endl;


    for(unsigned int i=1; i<Nt; i++)
    {
        //rhs=(1/dt*C.get_matrix()+lambda*Reaction)*solution.col(i-1)+C.get_rhs();
        //std::cout<<React.get_rhs()<<std::endl;
	rhs=(1/dt*C.get_matrix()-F_p.get_matrix()+F_m.get_matrix())*Ca.col(i-1)+F_p.get_rhs()-F_m.get_rhs();
        Ca.col(i)=solver.solve(rhs); 
    }
}



    
    

//Implicit transport upwind and esplicit reaction
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport,Data_2Reagents &data_2reagents) //As input there is the matrix solution where we store our solution at each istant,
//each row represent a spatial position, each column represent a time istant.
//In the vector vel there is the velocity evaluated at each node cell.
{ 
    //All the data that are needed to define the Transport System are extracted from the data structure
    auto &[L,phi,Nx,Nt,T,C_in,C_out,bc_cond]=data_transport;
    auto &[Ca_0,CaSiO3_0,K_sol,ph]=data_2reagents;

    //Computation of the spatial and temporal step from the data
    double h =static_cast<double>(L)/Nx;//constexpr?
    double dt=static_cast<double>(T)/Nt;//constexpr?

    //Porosità posta uguale a 1 
    //muparser_fun phi;
    //std::string s="1.0+0.0*x";
    //phi.set_value(s);

    //The Initial Condition are saved in an Eigen Vector. We recall that the value of the chemical species is saved in the middle of the cell (as the pressure in the Darcy System)
    for (unsigned int i=0; i<Nx; ++i)
        Ca(i,0)=Ca_0(h/2+i*h);
    //Ca(0,0)=1.0;

    Matrix_F_piu F_p(Nx,Nx);
    F_p.assemble_matrix(bc_cond,C_in,vel);

    Matrix_F_meno F_m(Nx,Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);

    Matrix_C C(Nx,Nx);
    C.assemble_matrix(phi,h);


    
    //Eigen::MatrixXd M(Nx,Nx);//Matrix of the linear system that has to be solved
    //Eigen::VectorXd rhs(Nx);//rhs of the lineay system that has to be solved

    //Here there is a different definition of the linear matrix of the system (and also of the rhs in the time loop)

    const Eigen::SparseMatrix<double> M(1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix());
   
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >   solver;
    // fill A and b;
    // Compute the ordering permutation vector from the structural pattern of A
    solver.analyzePattern(M); 
    // Compute the numerical factorization 
    solver.factorize(M); 
    
    Eigen::VectorXd rhs(Nx);

    for(unsigned int i=1; i<Nt; i++)
    {
        //rhs=(1/dt*C.get_matrix()+lambda*Reaction)*solution.col(i-1)+C.get_rhs();
        //std::cout<<React.get_rhs()<<std::endl;
	rhs=(1/dt*C.get_matrix())*Ca.col(i-1)+F_p.get_rhs()-F_m.get_rhs();
        Ca.col(i)=solver.solve(rhs);
    }
}



    
    
