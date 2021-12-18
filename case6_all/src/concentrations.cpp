#include "concentrations.hpp"
#include "matrix.hpp"
#include "gnuplot-iostream.hpp"
#include "functions.hpp"
#include <iostream>
#include <fstream>


Concentration::Concentration(const std::string &filename):data_transp(filename), data_reagents(filename), data_reaction(filename), data_CO2(filename)
{
    Ca=Matrix_full::Zero(data_transp.Nx,data_transp.Nt+1);
    H_piu=Matrix_full::Zero(data_transp.Nx,data_transp.Nt+1);
    HCO3_meno=Matrix_full::Zero(data_transp.Nx,data_transp.Nt+1);
    CO2=Matrix_full::Zero(data_transp.Nx,data_transp.Nt+1);
    CaSiO3=Matrix_full::Zero(data_transp.Nx,data_transp.Nt+1);
    SiO2=Matrix_full::Zero(data_transp.Nx,data_transp.Nt+1);

    if(data_reagents.method==1)
        method=ExplicitEuler;
    else if(data_reagents.method==2)
        method=PredictorCorrector;
    else
        method=Heun;
 
    h=static_cast<double>(data_transp.L/data_transp.Nx);
    dt=static_cast<double>(data_transp.T/data_transp.Nt);
}

unsigned int Concentration::get_Nx() const
{
    return data_transp.Nx;
}

unsigned int Concentration::get_Nt() const
{
    return data_transp.Nt;
}


void Concentration::set_initial_cond()
{

    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        Ca(i,0)=data_reagents.Ca_0(h+h*i);
        H_piu(i,0)=data_reagents.H_piu_0(h+h*i);
        HCO3_meno(i,0)=data_reagents.K_eq*data_reagents.CO2_0(h+h*i)/data_reagents.H_piu_0(h+h*i);
        CO2(i,0)=data_reagents.CO2_0(h+h*i);
        CaSiO3(i,0)=data_reagents.CaSiO3_0(h+h*i);
        SiO2(i,0)=data_reagents.SiO2_0(h+h*i);
    }
    
}

void Concentration::define_transport_solver(Solver& solver, Solver& solver1, Matrix& M_rhs, Vector& rhs_psi2, Vector& rhs_psi3, const Vector& vel, unsigned int Nx) 
{

    Matrix M(Nx,Nx);
    assemble_transport(M,M_rhs,vel);//function that assembles the transport part
    
    //Setting for the BC for the CO2 transport equation
    auto &[CO2_in,CO2_out,CO2_bc_cond,H_piu_in,H_piu_out,H_piu_bc_cond]=data_CO2;
    assemble_rhs(rhs_psi2,rhs_psi3,vel,CO2_in,CO2_out,CO2_bc_cond,H_piu_in,H_piu_out,H_piu_bc_cond);

    //Setting for the solid ODE equation    
    Matrix M_solid(Nx,Nx);//Mass matrix for solid (CaSiO3)
    Matrix rhs_solid(Nx,Nx);//rhs for the solid part
    assemble_transport(M_solid,rhs_solid,Vector::Zero(Nx+1));//function that assembles the mass matrix for the solid reagent
    
    set_solver(M,solver);
    
    set_solver(M_solid,solver1);
}


void Concentration::assemble_transport(Matrix& M, Matrix& M_rhs, const Vector& vel) const
{

    Matrix_C C(data_transp.Nx,data_transp.Nx);
    C.assemble_matrix(data_transp.phi,h);


    Matrix_F_piu F_p(data_transp.Nx,data_transp.Nx);
    F_p.set_data("In",0.0,vel);
    F_p.define_matrix();

    Matrix_F_meno F_m(data_transp.Nx,data_transp.Nx);
    F_m.set_data("In",0.0,vel);
    F_m.define_matrix();
 

    M=1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix();
 
    M_rhs=1/dt*C.get_matrix();

}


void Concentration::assemble_rhs(Vector& rhs_psi2, Vector& rhs_psi3, const Vector& vel, double CO2_in, double CO2_out, const std::string& CO2_bc, double H_piu_in, double H_piu_out, const std::string& H_piu_bc)
{
    Matrix_F_piu F_p_CO2(data_transp.Nx,data_transp.Nx);
    F_p_CO2.set_data(CO2_bc,CO2_in,vel);
    F_p_CO2.set_rhs();

    Matrix_F_meno F_m_CO2(data_transp.Nx,data_transp.Nx);
    F_m_CO2.set_data(CO2_bc,CO2_out,vel);
    F_m_CO2.set_rhs();
    
    if(CO2_bc=="In")
      CO2(0,0)=CO2_in;
    else
      CO2(data_transp.Nx-1,data_transp.Nt-1)=CO2_out;
      
    Matrix_F_piu F_p_H_piu(data_transp.Nx,data_transp.Nx);
    F_p_H_piu.set_data(H_piu_bc,H_piu_in,vel);
    F_p_H_piu.set_rhs();

    Matrix_F_meno F_m_H_piu(data_transp.Nx,data_transp.Nx);
    F_m_H_piu.set_data(H_piu_bc,H_piu_out,vel);
    F_m_H_piu.set_rhs();  
    
    if(H_piu_bc=="In")
      H_piu(0,0)=H_piu_in;
    else
      H_piu(data_transp.Nx-1,0)=H_piu_out;
    
    double C_in_HCO3_meno=data_reagents.K_eq*CO2_in/H_piu_in;
    double C_out_HCO3_meno=data_reagents.K_eq*CO2_out/H_piu_out;
   
    
    Matrix_F_piu F_p_HCO3_meno(data_transp.Nx,data_transp.Nx);
    F_p_HCO3_meno.set_data(CO2_bc,C_in_HCO3_meno,vel);
    F_p_HCO3_meno.set_rhs();

    Matrix_F_meno F_m_HCO3_meno(data_transp.Nx,data_transp.Nx);
    F_m_HCO3_meno.set_data(CO2_bc,C_out_HCO3_meno,vel);
    F_m_HCO3_meno.set_rhs();
    
    if(CO2_bc=="In")
      HCO3_meno(0,0)=C_in_HCO3_meno;
    else
      HCO3_meno(data_transp.Nx-1,0)=C_out_HCO3_meno;
    
    rhs_psi3=F_p_CO2.get_rhs()-F_m_CO2.get_rhs()+F_p_HCO3_meno.get_rhs()-F_m_HCO3_meno.get_rhs();

    rhs_psi2=F_m_HCO3_meno.get_rhs()-F_p_HCO3_meno.get_rhs()+F_p_H_piu.get_rhs()-F_m_H_piu.get_rhs();
    
}





void Concentration::compute_psi(unsigned int step, Vector& psi1, Vector& psi2, Vector& psi3, Vector& psi4, Vector& psi5) const
{
    psi1=Ca.col(step);
    psi2=H_piu.col(step)-HCO3_meno.col(step);
    psi3=CO2.col(step)+HCO3_meno.col(step);
    psi4=CaSiO3.col(step);
    psi5=SiO2.col(step);
}

void Concentration::compute_rd(unsigned int step, Vector& rd) const
{

    double temp;
    double omega=0.0;
    auto &[A, Rate_const, E, R, Temperature]=data_reaction;
    const double const_r= A*Rate_const*(std::exp(-E/(R*Temperature)));
    muparser_fun phi=data_transp.phi;

    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        temp=const_r*std::pow(H_piu(i,step),data_reagents.n);
        omega=Ca(i,step)*SiO2(i,step)/(H_piu(i,step)*H_piu(i,step));
        omega/=data_reagents.K_sol;
        rd(i)=phi(h/2+i*h)*temp*std::max((1-omega),0.);

    }
    
}

void Concentration::compute_rd_kd(unsigned int step, Vector& rd) const
{
    double omega;
    const double const_r=data_reagents.kd_i;
    muparser_fun phi=data_transp.phi;

    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        omega=Ca(i,step)*SiO2(i,step)/(H_piu(i,step)*H_piu(i,step));
        omega/=data_reagents.K_sol;
        rd(i)=phi(h/2+i*h)*const_r*std::max((1-omega),0.)*CaSiO3(i,step);
    }
}





void Concentration::one_step_transport_reaction(Vector& psi1, Vector& psi2, Vector& psi3, Vector& psi4, Vector& psi5, Vector& rd, const Matrix& M_rhs, const Vector& rhs_psi2, const Vector& rhs_psi3, unsigned int step, Solver &solver, Solver &solver1)
{
  switch(method)
    {
    case ExplicitEuler: //Explicit Euler
         Euler_Explicit(psi1,psi2,psi3,psi4,psi5,rd,M_rhs,rhs_psi2,rhs_psi3,solver,solver1);
      break;

    case PredictorCorrector: //Predictor Corrector
    {
        Vector psi1_{psi1}; //We need it because we need to store the temporary solution
        Vector psi2_{psi2};
        Vector psi3_{psi3};
        Vector psi4_{psi4};
        Vector psi5_{psi5};
     
        Euler_Explicit(psi1_,psi2_,psi3_,psi4_,psi5_,rd,M_rhs,rhs_psi2,rhs_psi3,solver,solver1);
        compute_concentration(step,psi1_,psi2_,psi3_,psi4_,psi5_);
        compute_rd_kd(step,rd);//new rd computation
  
        Euler_Explicit(psi1,psi2,psi3,psi4,psi5,rd,M_rhs,rhs_psi2,rhs_psi3,solver,solver1);
    }
    break;

    case Heun: //Heun
    {
        Vector psi1_{psi1};
        Vector psi2_{psi2};
        Vector psi3_{psi3};
        Vector psi4_{psi4};
        Vector psi5_{psi5};

        Euler_Explicit(psi1_,psi2_,psi3_,psi4_,psi5_,rd,M_rhs,rhs_psi2,rhs_psi3,solver,solver1);
        compute_concentration(step,psi1_,psi2_,psi3_,psi4_,psi5_);

        Vector rd_{Vector::Zero(data_transp.Nx)};
        compute_rd_kd(step,rd_);

        rd=0.5*rd+0.5*rd_;//new rd computation

        Euler_Explicit(psi1,psi2,psi3,psi4,psi5,rd,M_rhs,rhs_psi2,rhs_psi3,solver,solver1);
    }

    break;
    }

}

void Concentration::Euler_Explicit(Vector& psi1, Vector& psi2, Vector& psi3, Vector& psi4, Vector& psi5, const Vector& rd, const Matrix&  M_rhs, const Vector& rhs_psi2, const Vector& rhs_psi3, Solver &solver, Solver &solver1) const
{
    transport_and_reaction(psi1,M_rhs,Vector::Zero(psi1.size()),rd,solver);//reaction but not input bc (rhs=0) 
    transport_and_reaction(psi2,M_rhs,rhs_psi2,-2*rd,solver);//reaction with input bc
    transport_and_reaction(psi3,M_rhs,rhs_psi3,Vector::Zero(psi3.size()),solver);//not reaction but input boundary (rd=0)
    transport_and_reaction(psi4,M_rhs,Vector::Zero(psi4.size()),-rd,solver1);//reaction but not input bc (rhs=0) (different solver because solid is not moving)
    transport_and_reaction(psi5,M_rhs,Vector::Zero(psi5.size()),rd,solver);//reaction but not input bc (rhs=0)
}


void Concentration::transport_and_reaction(Vector& psi, const Matrix& M_rhs, const Vector& rhs, const Vector& rd, Solver &solver) const
{   
    const Vector temp{M_rhs*psi+rhs+rd*h};
    psi=solver.solve(temp);
}



void Concentration::compute_concentration(unsigned int step, const Vector& psi1, const Vector& psi2, const Vector& psi3, const Vector& psi4, const Vector& psi5)
{
    Vector old_it(6);//Initial guess::value at previous step
    Vector rhs(6);
    Matrix_full Jacob{Matrix_full::Zero(6,6)};
    Jacob(0,0)=1.0;
    Jacob(1,1)=1.0;
    Jacob(1,2)=-1.0;
    Jacob(2,2)=1.0;
    Jacob(2,3)=1.0;
    Jacob(3,4)=1.0;
    Jacob(4,5)=1.0;
    Jacob(5,3)=-data_reagents.K_eq;
    
    unsigned int max_iter=500;
    double tol=1.0e-15;

    for(unsigned int i=0; i<Ca.rows(); ++i)
    {
        old_it<<Ca(i,step-1),H_piu(i,step-1),HCO3_meno(i,step-1),CO2(i,step-1),CaSiO3(i,step-1),SiO2(i,step-1);  //initial guess
        
        double err=1;
        Vector dx(6);
        
        for(unsigned int iter=0; iter<max_iter and err>tol; ++iter)
        {
            compute_rhs(rhs, old_it, psi1(i), psi2(i), psi3(i), psi4(i), psi5(i));//F(x_k) computation;
            compute_Jacob(Jacob, old_it);//Jacobian computation
            const auto Jac=Jacob.fullPivLu();
            dx=Jac.solve(-rhs);//dx computation
            err=dx.norm()/old_it.norm();//error valutation
            old_it+=dx;//x_k+1           
        }
        Ca(i,step)=old_it(0);
        H_piu(i,step)=old_it(1);
        HCO3_meno(i,step)=old_it(2);
        CO2(i,step)=old_it(3);
        CaSiO3(i,step)=old_it(4);
        SiO2(i,step)=old_it(5);        
    }

}

void Concentration::compute_rhs(Vector& rhs, const Vector& old_it, double psi1, double psi2, double psi3, double psi4, double psi5) const
{

    double Ca{old_it(0)}, H_piu{old_it(1)}, HCO3_meno{old_it(2)}, CO2{old_it(3)}, CaSiO3{old_it(4)}, SiO2{old_it(5)};
    rhs<<Ca-psi1,H_piu-HCO3_meno-psi2,CO2+HCO3_meno-psi3,CaSiO3-psi4,SiO2-psi5,H_piu*HCO3_meno-data_reagents.K_eq*CO2;
}

void Concentration::compute_Jacob(Matrix_full& J,const Vector& old_it) const
{
    double HCO3_meno{old_it(2)},H_piu{old_it(1)};
    J(5,1)=HCO3_meno;
    J(5,2)=H_piu;
}


void Concentration::output_results_fixed_time(const std::string& name) const
{

    std::string filename;
    Matrix_full value1(data_transp.Nx,data_transp.Nt);
    if(name=="Ca")
    {
        filename="Ca_fixed_time.csv";
        value1=Ca;
    }
    else if(name=="H_piu")
    {
        filename="H_piu_fixed_time.csv";
        value1=H_piu;
    }
    else if(name=="HCO3_meno")
    {
        filename="HCO3_meno_fixed_time.csv";
        value1=HCO3_meno;
    }
    else if(name=="CO2")
    {
        filename="CO2_fixed_time.csv";
        value1=CO2;
    }
    else if(name=="CaSiO3")
    {
        filename="CaSiO3_fixed_time.csv";
        value1=CaSiO3;
    }
    else
    {
        filename="SiO2_fixed_time.csv";
        value1=SiO2;
    }

    std::ofstream file1(filename, std::ofstream::out);

    const unsigned int a=data_transp.Nt/10; //we print the result on the file each "a" seconds
      
    file1<<"space, ";
    for (unsigned int i=0.; i<data_transp.Nt+1; i+=a)
        file1<<"t="+std::to_string(i*dt)+"s, ";
    file1<<std::endl;

    const Vector x(Vector::LinSpaced(data_transp.Nx,h/2,data_transp.L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<data_transp.Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";

        for (unsigned int j=0; j<data_transp.Nt+1; j+=a)
        {
            file1<<value1(i,j)<<", ";
        }

        file1<<std::endl;

    }
    file1.close();

}


void Concentration::output_results_fixed_space(const std::string& name) const
{

    std::string filename;
    Matrix_full value1(data_transp.Nx,data_transp.Nt);
    if(name=="Ca")
    {
        filename="Ca_fixed_space.csv";
        value1=Ca;
    }
    else if(name=="H_piu")
    {
        filename="H_piu_fixed_space.csv";
        value1=H_piu;
    }
    else if(name=="HCO3_meno")
    {
        filename="HCO3_meno_fixed_space.csv";
        value1=HCO3_meno;
    }
    else if(name=="CO2")
    {
        filename="CO2_fixed_space.csv";
        value1=CO2;
    }
    else if(name=="CaSiO3")
    {
        filename="CaSiO3_fixed_space.csv";
        value1=CaSiO3;
    }
    else
    {
        filename="SiO2_fixed_space.csv";
        value1=SiO2;
    }
 
    std::ofstream file1(filename, std::ofstream::out);

    const unsigned int a=data_transp.Nx/10;//we print the result on the file each "a" meters

    file1<<"time, ";
    for (unsigned int i=0; i<data_transp.Nx; i+=a)
        file1<<"x="+std::to_string(h/2+i*h)+"m, ";
    file1<<"x="+std::to_string(data_transp.L-h/2)+"m, ";
    file1<<std::endl;

    const Vector t(Vector::LinSpaced(data_transp.Nt+1,0.0,data_transp.T));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<data_transp.Nt+1; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";

        for (unsigned int j=0; j<data_transp.Nx; j+=a)
        {
            file1<<value1(j,i)<<", ";

        }
        file1<<value1(data_transp.Nx-1,i)<<", ";
        file1<<std::endl;

    }
    file1.close();
}


void Concentration::output_all_reagents(unsigned int pos) const
{

    std::ofstream file1("all.csv", std::ofstream::out);
    file1<< "time, Ca, H+, HCO3-, CO2, CaSiO3, SiO2" << std::endl;

    const Vector t(Vector::LinSpaced(data_transp.Nt+1,0.0,data_transp.T));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<data_transp.Nt+1; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";

        file1<<Ca(pos,i)<<", ";
        file1<<H_piu(pos,i)<<", ";
        file1<<HCO3_meno(pos,i)<<", ";
        file1<<CO2(pos,i)<<", ";
        file1<<CaSiO3(pos,i)<<", ";
        file1<<SiO2(pos,i)<<", ";


        file1<<std::endl;

    }
    file1.close();

}

void Concentration::column(unsigned int i)
{
std::cout<<"Ca="<<Ca.col(i)<<std::endl;
std::cout<<"CaSiO3="<<CaSiO3.col(i)<<std::endl;
std::cout<<"CO2="<<CO2.col(i)<<std::endl;
std::cout<<"HCO3="<<HCO3_meno.col(i)<<std::endl;
std::cout<<"SiO2="<<SiO2.col(i)<<std::endl;
std::cout<<"H_piu="<<H_piu.col(i)<<std::endl;

}

