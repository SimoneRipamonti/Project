#include "concentrations.hpp"
#include "matrix.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>


Concentration::Concentration(const std::string &filename):data_transp(filename), data_reagents(filename), data_reaction(filename)
{
    Ca=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt);
    H_piu=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt);
    HCO3_meno=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt);
    CO2=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt);
    CaSiO3=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt);
    SiO2=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt);

    if(data_reagents.method==1)
        method=EsplicitEuler;
    else if(data_reagents.method==2)
        method=PredictorCorrector;
    else
        method=Heun;
}

unsigned int Concentration::get_Nx()
{
    return data_transp.Nx;
}

unsigned int Concentration::get_Nt()
{
    return data_transp.Nt;
}


void Concentration::set_initial_cond()
{

    double h{static_cast<double>(data_transp.L/data_transp.Nx)};
    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        Ca(i,0)=data_reagents.Ca_0(h+h*i);
        H_piu(i,0)=data_reagents.H_piu_0(h+h*i);
        HCO3_meno(i,0)=4.45e-7*data_reagents.CO2_0(h+h*i)/data_reagents.H_piu_0(h+h*i);
        CO2(i,0)=data_reagents.CO2_0(h+h*i);
        CaSiO3(i,0)=data_reagents.CaSiO3_0(h+h*i);
        SiO2(i,0)=data_reagents.SiO2_0(h+h*i);
    }


}

void Concentration::assemble_transport(Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& vel)
{

    double h{static_cast<double>(data_transp.L/data_transp.Nx)};
    double dt=static_cast<double>(data_transp.T/data_transp.Nt);


    Matrix_C C(data_transp.Nx,data_transp.Nx);
    C.assemble_matrix(data_transp.phi,h);


    Matrix_F_piu F_p(data_transp.Nx,data_transp.Nx);

    F_p.assemble_matrix("In",0.0,vel);
    Matrix_F_meno F_m(data_transp.Nx,data_transp.Nx);
    F_m.assemble_matrix("In",0.0,vel);


    M=1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix();


    rhs=1/dt*C.get_matrix();

}




void Concentration::assemble_transport_CO2(Eigen::SparseMatrix<double>& M_CO2, Eigen::VectorXd& rhs_CO2, const Eigen::VectorXd& vel, double C_in, double C_out, const std::string& bc_cond)
{

    double h{static_cast<double>(data_transp.L/data_transp.Nx)};
    double dt=static_cast<double>(data_transp.T/data_transp.Nt);

    Matrix_C C(data_transp.Nx,data_transp.Nx);
    C.assemble_matrix(data_transp.phi,h);


    Matrix_F_piu F_p(data_transp.Nx,data_transp.Nx);

    F_p.assemble_matrix(bc_cond,C_in,vel);
    Matrix_F_meno F_m(data_transp.Nx,data_transp.Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);



    M_CO2=1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix();

    rhs_CO2=F_p.get_rhs()-F_m.get_rhs();

}




void Concentration::compute_phi(unsigned int step, Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, Eigen::VectorXd& phi5)
{
    phi1=Ca.col(step);
    phi2=H_piu.col(step)-HCO3_meno.col(step);
    phi3=CO2.col(step)+HCO3_meno.col(step);
    phi4=CaSiO3.col(step);
    phi5=SiO2.col(step);
}

void Concentration::compute_rd(unsigned int step, Eigen::VectorXd& rd)
{

    //double temp;
    double omega;
    auto [A, Rate_const, E, R, Temperature]=data_reaction;
    const double const_r= A*Rate_const*(std::exp(-E/(R*Temperature)));

    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        //temp=const_r*std::pow(H_piu(i,step),data_reagents.n);

        omega=Ca(i,step)*SiO2(i,step)/(H_piu(i,step)*H_piu(i,step));//nel codice di Anna c'è anche un Hpiu al quadrato
        omega/=data_reagents.K_sol;
        //rd(i)=temp*std::max((1-omega),0.);
        rd(i)=const_r*std::max((1-omega),0.);

    }

}

void Concentration::compute_rd_kp(unsigned int step, Eigen::VectorXd& rd)
{

    double omega;
    const double const_r= data_reaction.A*data_reagents.kp_i;

    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        omega=Ca(i,step)*SiO2(i,step)/(H_piu(i,step)*H_piu(i,step));//nel codice di Anna c'è anche un Hpiu al quadrato
        omega/=data_reagents.K_sol;
        rd(i)=const_r*std::max((1-omega),0.);
    }

}





void Concentration::one_step_transport_reaction(Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5, Eigen::VectorXd& rd, const  Eigen::SparseMatrix<double> rhs, const Eigen::VectorXd& rhs_CO2, unsigned int step, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2)
{

    switch(method)
    {
    case EsplicitEuler: //Euler Esplicit

        //Euler_Esplicit(phi1,phi2,phi3,phi4,phi5,rd,M,rhs,solver);
        Euler_Esplicit(psi1,psi2,psi3,psi4,psi5,rd,rhs,rhs_CO2,solver,solver2);
        break;

    case PredictorCorrector: //Predictor Corrector
    {
        Eigen::VectorXd psi1_{psi1}; //Mi servono perché devo conservare le phi al passo n
        Eigen::VectorXd psi2_{psi2};
        Eigen::VectorXd psi3_{psi3};
        Eigen::VectorXd psi4_{psi4};
        Eigen::VectorXd psi5_{psi5};

        Euler_Esplicit(psi1_,psi2_,psi3_,psi4_,psi5_,rd,rhs,rhs_CO2,solver,solver2);
        compute_concentration(step,psi1_,psi2_,psi3_,psi4_,psi5_);
        compute_rd(step,rd);//calcolo il nuovo rd


        //P_C(phi1,phi2,phi3,phi4,phi5,rd,rhs,solver);
        Euler_Esplicit(psi1,psi2,psi3,psi4,psi5,rd,rhs,rhs_CO2,solver,solver2);
    }
    break;

    case Heun: //Heun
    {
        Eigen::VectorXd psi1_{psi1};
        Eigen::VectorXd psi2_{psi2};
        Eigen::VectorXd psi3_{psi3};
        Eigen::VectorXd psi4_{psi4};
        Eigen::VectorXd psi5_{psi5};

        Euler_Esplicit(psi1_,psi2_,psi3_,psi4_,psi5_,rd,rhs,rhs_CO2,solver,solver2);
        compute_concentration(step,psi1_,psi2_,psi3_,psi4_,psi5_);

        Eigen::VectorXd rd_{Eigen::VectorXd::Zero(data_transp.Nx)};
        compute_rd(step,rd_);

        rd=0.5*rd+0.5*rd_;//ottengo il nuovo rd complessivo

        Euler_Esplicit(psi1,psi2,psi3,psi4,psi5,rd,rhs,rhs_CO2,solver,solver2);
    }

    break;
    }

}

void Concentration::Euler_Esplicit(Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5, const Eigen::VectorXd& rd, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rhs_CO2, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2)
{
    transport_and_reaction(psi1,rhs,rd,solver);
    transport_and_reaction(psi2,rhs,Eigen::VectorXd::Zero(psi2.size()),solver);//no reaction
    //transport_and_reaction(phi2,rhs,-2*rd,solver);
    transport_and_reaction_CO2(psi3,rhs,rhs_CO2,Eigen::VectorXd::Zero(psi2.size()),solver2);//no reaction
    transport_and_reaction(psi4,rhs,-rd,solver);
    //transport_and_reaction(phi5,M,rhs,rd,solver);
    transport_and_reaction(psi5,rhs,rd,solver);
}



void Concentration::transport_and_reaction(Eigen::VectorXd& psi, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rd,Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver)
{

    double h{static_cast<double>(data_transp.L/data_transp.Nx)};
    const Eigen::VectorXd temp=rhs*psi+rd*h;
    psi=solver.solve(temp);
}

void Concentration::transport_and_reaction_CO2(Eigen::VectorXd& psi, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rhs_CO2, const Eigen::VectorXd& rd, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2)
{
    double h{static_cast<double>(data_transp.L/data_transp.Nx)};
    const Eigen::VectorXd temp=rhs*psi+rhs_CO2+rd*h;
    psi=solver2.solve(temp);
}



void Concentration::compute_concentration(unsigned int step, const Eigen::VectorXd& psi1, const Eigen::VectorXd& psi2, const Eigen::VectorXd& psi3, const Eigen::VectorXd& psi4, const Eigen::VectorXd& psi5)
{

//RISOLVO IL SISTEMA NON LINEARE CON NEWTON
    //Initial guess::value at previous step
    Eigen::VectorXd old_it(6);
    Eigen::VectorXd rhs(6);
    Eigen::MatrixXd Jacob{Eigen::MatrixXd::Zero(6,6)};
    Jacob(0,0)=1.0;
    Jacob(1,1)=1.0;
    Jacob(1,2)=-1.0;
    Jacob(2,2)=1.0;
    Jacob(2,3)=1.0;
    Jacob(3,4)=1.0;
    Jacob(4,5)=1.0;
    Jacob(5,3)=-data_reagents.K_eq;

    for(unsigned int i=0; i<Ca.rows(); ++i)
    {
        old_it<<Ca(i,step-1),H_piu(i,step-1),HCO3_meno(i,step-1),CO2(i,step-1),CaSiO3(i,step-1),SiO2(i,step-1);  //guess iniziale

        unsigned int max_iter=500;
        double tol=1.0e-14;
        double err=1;
        Eigen::VectorXd dx(6);

        for(unsigned int iter=0; iter<max_iter and err>tol; ++iter)
        {
            compute_rhs(rhs, old_it, psi1(i), psi2(i), psi3(i), psi4(i), psi5(i));///Calcolo F(x_k);
            compute_Jacob(Jacob, old_it);//Calcolo Jacob
            const auto Jac=Jacob.fullPivLu();
            dx=Jac.solve(-rhs);//Calcolo dx;
            err=dx.norm();//Valuto errore
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

void Concentration::compute_rhs(Eigen::VectorXd& rhs,const Eigen::VectorXd& old_it, double psi1, double psi2, double psi3, double psi4, double psi5)
{

    double Ca{old_it(0)}, H_piu{old_it(1)}, HCO3_meno{old_it(2)}, CO2{old_it(3)}, CaSiO3{old_it(4)}, SiO2{old_it(5)};
    rhs<<Ca-psi1,H_piu-HCO3_meno-psi2,CO2+HCO3_meno-psi3,CaSiO3-psi4,SiO2-psi5,H_piu*HCO3_meno-data_reagents.K_eq*CO2;
}

void Concentration::compute_Jacob(Eigen::MatrixXd& J,const Eigen::VectorXd& old_it)
{
    double HCO3_meno{old_it(2)},H_piu{old_it(1)};
    J(5,1)=HCO3_meno;
    J(5,2)=H_piu;
}


void Concentration::output_results_fixed_time(const std::string& name)
{

    std::string filename;
    Eigen::MatrixXd value1(data_transp.Nx,data_transp.Nt);
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

    file1<< "space, t0,...,t_Nt-1" << std::endl;

    double h =static_cast<double>(data_transp.L)/data_transp.Nx;
    const Eigen::VectorXd x(Eigen::VectorXd::LinSpaced(data_transp.Nx,h/2,data_transp.L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<data_transp.Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";

        for (unsigned int j=0; j<data_transp.Nt; ++j)
        {
            file1<<value1(i,j)<<", ";
        }

        file1<<std::endl;

    }
    file1.close();

}


void Concentration::output_results_fixed_space(const std::string& name)
{

    std::string filename;
    Eigen::MatrixXd value1(data_transp.Nx,data_transp.Nt);
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
    file1<< "time, x0,...,x_Nx-1" << std::endl;


    double dt=static_cast<double>(data_transp.T)/data_transp.Nt;
    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(data_transp.Nt,0.0,data_transp.T-dt));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<data_transp.Nt; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< t[i] <<", ";

        for (unsigned int j=0; j<data_transp.Nx; ++j)
        {
            file1<<value1(j,i)<<", ";

        }
        file1<<std::endl;

    }
    file1.close();
}


void Concentration::output_all_reagents(unsigned int pos)
{

    std::ofstream file1("all.csv", std::ofstream::out);
    file1<< "time, Ca, H+, HCO3-, CO2, CaSiO3, SiO2" << std::endl;


    double dt=static_cast<double>(data_transp.T)/data_transp.Nt;
    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(data_transp.Nt,0.0,data_transp.T-dt));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<data_transp.Nt; ++i) //Loop to save the matrix by column in the CSV file
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




