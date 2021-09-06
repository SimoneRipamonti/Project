#include "concentrations.hpp"
#include "matrix.hpp"
#include "gnuplot-iostream.hpp"
#include <iostream>
#include <fstream>


Concentration::Concentration(const std::string &filename):data_transp(filename), data_reagents(filename), data_reaction(filename)
{
    Ca=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt+1);
    H_piu=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt+1);
    HCO3_meno=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt+1);
    CO2=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt+1);
    CaSiO3=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt+1);
    SiO2=Eigen::MatrixXd::Zero(data_transp.Nx,data_transp.Nt+1);

    if(data_reagents.method==1)
        method=EsplicitEuler;
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
        //HCO3_meno(i,0)=data_reagents.K_eq*data_reagents.CO2_0(h+h*i)/data_reagents.H_piu_0(h+h*i);
        HCO3_meno(i,0)=0.0;
        CO2(i,0)=data_reagents.CO2_0(h+h*i);
        CaSiO3(i,0)=data_reagents.CaSiO3_0(h+h*i);
        SiO2(i,0)=data_reagents.SiO2_0(h+h*i);
    }


}

void Concentration::assemble_transport(Eigen::SparseMatrix<double>& M, Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& vel) const
{


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

    Matrix_C C(data_transp.Nx,data_transp.Nx);
    C.assemble_matrix(data_transp.phi,h);


    Matrix_F_piu F_p(data_transp.Nx,data_transp.Nx);

    F_p.assemble_matrix(bc_cond,C_in,vel);
    Matrix_F_meno F_m(data_transp.Nx,data_transp.Nx);
    F_m.assemble_matrix(bc_cond,C_out,vel);
    
    if(bc_cond=="In")
      CO2(0,0)=C_in;
    else
      CO2(data_transp.Nx-1,data_transp.Nt-1)=C_out;
    

   //std::cout<<"F_p="<<std::endl;
   //std::cout<<F_p.get_matrix()<<std::endl;
   //std::cout<<"F_m="<<std::endl;
   //std::cout<<F_m.get_matrix()<<std::endl;



    M_CO2=1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix();

    //std::cout<<M_CO2<<std::endl;

    rhs_CO2=F_p.get_rhs()-F_m.get_rhs();

}




void Concentration::compute_psi(unsigned int step, Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5) const
{
   /* psi1=Ca.col(step);
    //psi2=H_piu.col(step)-HCO3_meno.col(step);
    psi2=H_piu.col(step);
    psi3=CO2.col(step);
    psi4=CaSiO3.col(step);
    psi5=SiO2.col(step);*/

    
    psi1=Ca.col(step);
    psi2=H_piu.col(step)-HCO3_meno.col(step);
    psi3=CO2.col(step)+HCO3_meno.col(step);
    psi4=CaSiO3.col(step);
    psi5=SiO2.col(step);
}

void Concentration::compute_rd(unsigned int step, Eigen::VectorXd& rd) const
{

    double temp;
    double omega;
    auto [A, Rate_const, E, R, Temperature]=data_reaction;
    const double const_r= A*Rate_const*(std::exp(-E/(R*Temperature)));
    /*if(step==0)
      {std::cout<<"const_r="<<const_r<<std::endl;}*/
    muparser_fun phi=data_transp.phi;

    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        temp=const_r*std::pow(H_piu(i,step),data_reagents.n);

        //omega=Ca(i,step)*SiO2(i,step)/(H_piu(i,step)*H_piu(i,step));//nel codice di Anna c'è anche un Hpiu al quadrato
        omega/=data_reagents.K_sol;
        /*if(step==0 && i==0)
        {std::cout<<"omega="<<omega<<std::endl;
         std::cout<<"temp="<<temp<<std::endl;}*/

        rd(i)=phi(h/2+i*h)*temp*std::max((1-omega),0.)*CaSiO3(i,step);
        //rd(i)=phi(h/2+i*h)*temp*std::max((1-omega),0.);

    }
    
    /*if(step==0)
    {unsigned int a=data_transp.Nx/10;
    std::cout<<"rd="<<std::endl;
       for(unsigned int j=0;j<data_transp.Nx;j+=a)
           {std::cout<<rd(j)<<std::endl;}
    }*/

}

void Concentration::compute_rd_kp(unsigned int step, Eigen::VectorXd& rd) const
{

    double omega;
    const double const_r= data_reaction.A*data_reagents.kp_i;
    muparser_fun phi=data_transp.phi;

    for (unsigned int i=0; i<Ca.rows(); ++i)
    {
        omega=Ca(i,step)*SiO2(i,step)/(H_piu(i,step)*H_piu(i,step));//nel codice di Anna c'è anche un Hpiu al quadrato
        omega/=data_reagents.K_sol;
        rd(i)=phi(h/2+i*h)*const_r*std::max((1-omega),0.);
    }


}





void Concentration::one_step_transport_reaction(Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5, Eigen::VectorXd& rd, const  Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rhs_CO2, unsigned int step, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2) 
{

 /*if(step==1)
    {  
       //std::cout<<rhs<<std::endl;
       const Eigen::VectorXd temp1=rhs*psi1;//+rd*h;
       unsigned int a=data_transp.Nx/10;
       std::cout<<"temp1="<<std::endl;
       for(unsigned int j=0;j<data_transp.Nx;j+=a)
             {std::cout<<temp1(j)<<std::endl;}

        const Eigen::VectorXd temp2=rhs*psi2;
       
       std::cout<<"temp2="<<std::endl;
       for(unsigned int j=0;j<data_transp.Nx;j+=a)
             {std::cout<<temp2(j)<<std::endl;}

         const Eigen::VectorXd temp3=rhs*psi3;
       
       std::cout<<"temp3="<<std::endl;
       for(unsigned int j=0;j<data_transp.Nx;j+=a)
             {std::cout<<temp3(j)<<std::endl;}

         const Eigen::VectorXd temp4=rhs*psi4-rd*h;
       
       std::cout<<"temp4="<<std::endl;
       for(unsigned int j=0;j<data_transp.Nx;j+=a)
             {std::cout<<temp4(j)<<std::endl;}

         const Eigen::VectorXd temp5=rhs*psi5+rd*h;
       
       std::cout<<"temp5="<<std::endl;
       for(unsigned int j=0;j<data_transp.Nx;j+=a)
             {std::cout<<temp5(j)<<std::endl;}

     


    }*/




    
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

void Concentration::Euler_Esplicit(Eigen::VectorXd& psi1, Eigen::VectorXd& psi2, Eigen::VectorXd& psi3, Eigen::VectorXd& psi4, Eigen::VectorXd& psi5, const Eigen::VectorXd& rd, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rhs_CO2, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2) const
{
    transport_and_reaction(psi1,rhs,rd,solver);
    //transport_and_reaction(psi2,rhs,Eigen::VectorXd::Zero(psi2.size()),solver);//no reaction
    transport_and_reaction(psi2,rhs,-2*rd,solver);
    transport_and_reaction_CO2(psi3,rhs,rhs_CO2,Eigen::VectorXd::Zero(psi2.size()),solver2);//no reaction
    transport_and_reaction(psi4,rhs,-rd,solver);
    //transport_and_reaction(phi5,M,rhs,rd,solver);
    transport_and_reaction(psi5,rhs,rd,solver);
}



void Concentration::transport_and_reaction(Eigen::VectorXd& psi, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rd,Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver) const
{
    const Eigen::VectorXd temp=rhs*psi+rd*h;
    psi=solver.solve(temp);
}

void Concentration::transport_and_reaction_CO2(Eigen::VectorXd& psi, const Eigen::SparseMatrix<double>& rhs, const Eigen::VectorXd& rhs_CO2, const Eigen::VectorXd& rd, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver2) const
{   
    const Eigen::VectorXd temp=rhs*psi+rhs_CO2+rd*h;
    psi=solver2.solve(temp);
}



void Concentration::compute_concentration(unsigned int step, const Eigen::VectorXd& psi1, const Eigen::VectorXd& psi2, const Eigen::VectorXd& psi3, const Eigen::VectorXd& psi4, const Eigen::VectorXd& psi5)
{

/*for(unsigned int i=0; i<Ca.rows(); ++i)

{       Ca(i,step)=psi1(i);
        H_piu(i,step)=psi2(i);
        CO2(i,step)=psi3(i);
        CaSiO3(i,step)=psi4(i);
        SiO2(i,step)=psi5(i);

}*/


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

void Concentration::compute_rhs(Eigen::VectorXd& rhs,const Eigen::VectorXd& old_it, double psi1, double psi2, double psi3, double psi4, double psi5) const
{

    double Ca{old_it(0)}, H_piu{old_it(1)}, HCO3_meno{old_it(2)}, CO2{old_it(3)}, CaSiO3{old_it(4)}, SiO2{old_it(5)};
    rhs<<Ca-psi1,H_piu-HCO3_meno-psi2,CO2+HCO3_meno-psi3,CaSiO3-psi4,SiO2-psi5,H_piu*HCO3_meno-data_reagents.K_eq*CO2;
}

void Concentration::compute_Jacob(Eigen::MatrixXd& J,const Eigen::VectorXd& old_it) const
{
    double HCO3_meno{old_it(2)},H_piu{old_it(1)};
    J(5,1)=HCO3_meno;
    J(5,2)=H_piu;
}


void Concentration::output_results_fixed_time(const std::string& name) const
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

    const unsigned int a=data_transp.Nt/10; //we print the result on the file each "a" seconds

    file1<<"space, ";
    for (unsigned int i=0.; i<data_transp.Nt+1; i+=a)
        file1<<"t="+std::to_string(i*dt)+"s, ";
    file1<<std::endl;

    const Eigen::VectorXd x(Eigen::VectorXd::LinSpaced(data_transp.Nx,h/2,data_transp.L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

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

    const unsigned int a=data_transp.Nx/10;//we print the result on the file each "a" meters

    file1<<"time, ";
    for (unsigned int i=0; i<data_transp.Nx; i+=a)
        file1<<"x="+std::to_string(h/2+i*h)+"m, ";
    file1<<"x="+std::to_string(data_transp.L-h/2)+"m, ";
    file1<<std::endl;

    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(data_transp.Nt+1,0.0,data_transp.T));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

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

    const Eigen::VectorXd t(Eigen::VectorXd::LinSpaced(data_transp.Nt+1,0.0,data_transp.T));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

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

/*
void Concentration::output_results_fixed_time_psi(const std::string& name, const Eigen::VectorXd& psi, unsigned step) const
{

    std::string filename;
    //Eigen::MatrixXd value1(data_transp.Nx,data_transp.Nt);
    if(name=="psi1")
    {
        filename="psi1_fixed_time.csv";
        //value1=Ca;
    }
    else if(name=="psi2")
    {
        filename="psi2_fixed_time.csv";
        //value1=H_piu;
    }
    else if(name=="psi3")
    {
        filename="psi3_fixed_time.csv";
        //value1=HCO3_meno;
    }
    else if(name=="CO2")
    {
        filename="psi4_fixed_time.csv";
        //value1=CO2;
    }
    else
    {
        filename="psi5_fixed_time.csv";
        //value1=SiO2;
    }

    std::ofstream file1(filename, std::ofstream::out);

    const unsigned int a=data_transp.Nt/10; //we print the result on the file each "a" seconds
    file1<<"space, ";
    for (unsigned int i=0.; i<data_transp.Nt; i+=a)
        file1<<"t="+std::to_string(i*dt)+"s, ";
    file1<<std::endl;

    const Eigen::VectorXd x(Eigen::VectorXd::LinSpaced(data_transp.Nx,h/2,data_transp.L-h/2));//Definition of the space vector (Concnetration values are stored in the middle of the cell)

    for (unsigned int i = 0; i<data_transp.Nx; ++i) //Loop to save the matrix by column in the CSV file
    {
        file1<< x[i] <<", ";

        file1<<value1(i,j)<<", ";
        

        file1<<std::endl;

    }
    file1.close();

}*/





