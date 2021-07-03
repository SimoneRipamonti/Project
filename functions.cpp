#include "functions.hpp"

void set_initial_cond(Eigen::VectorXd& Ca, Eigen::VectorXd& H_piu, Eigen::VectorXd& HCO3_meno, Eigen::VectorXd& CO2, Eigen::VectorXd& CaSiO3, Eigen::VectorXd& SiO2, muparser_fun Ca0, muparser_fun H_piu0,
                      muparser_fun HCO3_meno0, muparser_fun CO20, muparser_fun CaSiO30, muparser_fun SiO20, double h)
{
 for (unsigned int i=0; Ca.size();++i)
  {
   Ca(i)=Ca0(h+h*i);
   H_piu(i)=H_piu0(h+h*i);
   HCO3_meno(i)=HCO3_meno0(h+h*i);
   CO2(i)=CO20(h+h*i);
   CaSiO3(i)=CaSiO30(h+h*i);
   SiO20(i)=SiO20(h+h*i);
  }
}


void assemble_transport(Eigen::MatrixXd& M, Eigen::VectorXd& rhs, Eigen::VectorXd vel, muparser_fun phi, double h)
{
   Matrix_F_piu F_p(Nx,Nx);	
   F_p.assemble_matrix("In",0.0,vel);
 
   Matrix_F_meno F_m(Nx,Nx);
   F_m.assemble_matrix("In",0.0,vel);

   Matrix_C C(Nx,Nx);
   C.assemble_matrix(phi,h);

   M=1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix();
	
   rhs=(1/dt*C.get_matrix());
}

//calcolo phi
void compute_phi(Eigen::VectorXd& phi1,Eigen::VectorXd& phi2,Eigen::VectorXd& phi3,Eigen::VectorXd& phi4,Eigen::VectorXd& phi5,Eigen::VectorXd Ca,Eigen::VectorXd H_piu,Eigen::VectorXd HCO3_meno, Eigen::VectorXd CO2,Eigen::VectorXd CaSiO3,Eigen::VectorXd SiO2)
{
  phi1=Ca;
  phi2=H_piu-HCO3_meno;
  phi3=CO2+HCO3_meno;
  phi4=CaSiO3;
  phi5=SiO2;
}

//calcolo rd
void compute_rd(Eigen::VectorXd& rd, const Eigen::VectorXd Ca, const Eigen::VectorXd H_piu, const Eigen::VectorXd SiO2,const double const_r, double K_eq, double n) 
{
  const Eigen::VectorXd temp=const_r*H_piu.array().pow(n);
  const Eigen::VectorXd omega=(Ca.cwiseProduct(SiO2)).cwiseQuotient(H_piu)/K_eq;
  rd=temp.cwiseProduct(Eigen::VectorXd::Ones(Ca.size())-omega); 
}

//one_step_reaction
void one_step_transport_reaction(Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, Eigen::VectorXd& phi5, const Eigen::VectorXd rd, const auto M_lu, Eigen::MatrixXd rhs)
{
  transport_and_reaction(&phi1,M_lu,rhs,rd);
  transport_and_reaction(&phi2,M_lu,rhs,Eigen::Vector::Zeros(phi2.size());//no reaction
  transport_and_reaction(&phi3,M_lu,rhs,Eigen::Vector::Zeros(phi3.size());//no reaction
  transport_and_reaction(&phi4,M_lu,rhs,-rd);
  transport_and_reaction(&phi5,M_lu,rhs,rd);
}

//Trasporto e reazione
void transport_and_reaction(Eigen::VectorXd& phi, const auto M_lu, const Eigen::MatrixXd &rhs, const Eigen::VectorXd rd)
{
  const Eigen::VectorXd temp=rhs*phi+rd;
  phi=M_lu.solve(temp);
}


//Risoluzione sistema non lineare
compute_concentration(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, Eigen::MatrixXd& SiO2, const Eigen::VectorXd phi1, const Eigen::VectorXd phi2, const Eigen::VectorXd phi3, const Eigen::VectorXd phi4, const Eigen::VectorXd phi5, unsigned int step, double K_eq)
{

  //RISOLVO IL SISTEMA NON LINEARE CON NEWTON
  //Initial guess::value at previous step
  Eigen::VectorXd old_it(6);
  Eigen::VectorXd rhs(6);
  Eigen::MatrixXd::Zeros(6,6) Jacob;
  Jacob(0,0)=1.0;
  Jacob(1,1)=1.0;
  Jacob(1,2)=-1.0;
  Jacob(2,2)=1.0;
  Jacob(2,3)=1.0;
  Jacob(3,4)=1.0;
  Jacob(4,5)=1.0;
  Jacob(5,3)=-K_eq;

  for(unsigned int i=0; i<Nx; ++i)
  {
   old_it<<Ca(i,step-1),H_piu(i,step-1),HCO3_meno(i,step-1),CO2(i,step-1),CaSiO3(i,step-1),SiO2(i,step-1); //guess iniziale
   
   unsigned int max_iter=500;
   double tol=1.0e-14;
   double err=1; 
   Eigen::VectorXd dx(6);

      for(unsigned int iter=0; iter<max_iter and err>tol; ++iter)
       {
        compute_rhs(rhs, old_it, phi1(i), phi2(i), phi3(i), phi4(i), phi5(i),K_eq);///Calcolo F(x_k);
        compute_Jacob_last_row(Jacob, old_it);//Calcolo Jacob
        dx=Jacob.solve(-rhs);//Calcolo dx;
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
  

void compute_rhs(Eigen::VectorXd& rhs,const Eigen::VectorXd old_it, double phi1, double phi2, double phi3, double phi4, double phi5, double K_eq)
{

  double Ca(old_it(0)), H_piu(old_it(1)), HCO3_meno(old_it(2)), CO2(old_it(3)), CaSiO3(old_it(4)),   SiO2(old_it(5));
  rhs<<Ca-phi1,H_piu-HCO3_meno-phi2,CO2+HCO3_meno-phi3,CaSiO3-phi4,SiO2-phi5,H_piu*HCO3_meno/CO2-K_eq;

}


void compute_Jacob(Eigen::MatrixXd& J,const Eigen::VectorXd old_it)
{
     double HCO3_meno(old_it(2)),H_piu(old_it(1));
     J(5,1)=HCO3_meno;
     J(5,2)=H_piu;
}

