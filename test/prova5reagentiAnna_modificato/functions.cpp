#include "functions.hpp"


void change_element (Eigen::MatrixXd& a){
a(0,0)=1;
}

void set_initial_cond(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, muparser_fun Ca0, muparser_fun H_piu0,
                      muparser_fun HCO3_meno0, muparser_fun CO20, muparser_fun CaSiO30, double h)
{
 for (unsigned int i=0; i<Ca.rows();++i)
  {
   Ca(i,0)=Ca0(h+h*i);
   H_piu(i,0)=H_piu0(h+h*i);
   HCO3_meno(i,0)=4.45e-7*CO20(h+h*i)/H_piu0(h+h*i);
   CO2(i,0)=CO20(h+h*i);
   CaSiO3(i,0)=CaSiO30(h+h*i);
   //SiO20(i,0)=SiO20(h+h*i);
  }
}



void assemble_transport(Eigen::SparseMatrix<double>& M, Eigen::MatrixXd& rhs, Eigen::VectorXd vel, muparser_fun phi, double h, unsigned int Nx, double dt)
{
   Matrix_C C(Nx,Nx);
   C.assemble_matrix(phi,h);

   
   Matrix_F_piu F_p(Nx,Nx);
   
   F_p.assemble_matrix("In",0.0,vel);
   Matrix_F_meno F_m(Nx,Nx);
   F_m.assemble_matrix("In",0.0,vel);

   
   M=(1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix()).sparseView();

	
   rhs=1/dt*C.get_matrix();
 
}

//calcolo phi
void compute_phi(Eigen::VectorXd& phi1,Eigen::VectorXd& phi2,Eigen::VectorXd& phi3,Eigen::VectorXd& phi4, const Eigen::VectorXd Ca, const Eigen::VectorXd H_piu, const Eigen::VectorXd HCO3_meno, const Eigen::VectorXd CO2, const Eigen::VectorXd CaSiO3)
{
  phi1=Ca;
  phi2=H_piu-HCO3_meno;
  phi3=CO2+HCO3_meno;
  phi4=CaSiO3;
  //phi5=SiO2;
}

//calcolo rd
void compute_rd(Eigen::VectorXd& rd, const Eigen::VectorXd Ca, const Eigen::VectorXd H_piu, double const_r, double K_sol, double n) 
{
  /*const Eigen::VectorXd temp=const_r*H_piu.array().pow(n);
  const Eigen::VectorXd omega=(Ca.cwiseProduct(SiO2)).cwiseQuotient(H_piu)/K_eq;
  rd=temp.cwiseProduct(Eigen::VectorXd::Ones(Ca.size())-omega);*/ 
  ///
  double temp;
  double omega;
  for (unsigned int i=0; i<Ca.size(); ++i)
   {   
       //temp=const_r*std::pow(H_piu(i),n);
       temp=const_r;
       omega=Ca(i)*Ca(i)/(H_piu(i)*H_piu(i));
       omega/=K_sol;
       rd(i)=temp*std::max((1-omega),0.);
    }

}

//one_step_reaction
void one_step_transport_reaction(Eigen::VectorXd& phi1, Eigen::VectorXd& phi2, Eigen::VectorXd& phi3, Eigen::VectorXd& phi4, const Eigen::VectorXd rd, const Eigen::MatrixXd rhs, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver)
{
  transport_and_reaction(phi1,rhs,rd,solver);
  transport_and_reaction(phi2,rhs,-2*rd,solver);//no reaction
  transport_and_reaction(phi3,rhs,Eigen::VectorXd::Zero(phi3.size()),solver);//no reaction
  transport_and_reaction(phi4,rhs,-rd,solver);
  //transport_and_reaction(phi5,M_lu,rhs,rd);
}

//Trasporto e reazione
void transport_and_reaction(Eigen::VectorXd& phi, const Eigen::MatrixXd rhs, const Eigen::VectorXd rd, Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > &solver)
{ 
 const Eigen::VectorXd temp=rhs*phi+rd;
 phi=solver.solve(temp);
}



//Risoluzione sistema non lineare
void compute_concentration(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, const Eigen::VectorXd phi1, const Eigen::VectorXd phi2, const Eigen::VectorXd phi3, const Eigen::VectorXd phi4, unsigned int step, double K_eq)
{

  //RISOLVO IL SISTEMA NON LINEARE CON NEWTON
  //Initial guess::value at previous step
  Eigen::VectorXd old_it(5);
  Eigen::VectorXd rhs(5);
  Eigen::MatrixXd Jacob{Eigen::MatrixXd::Zero(5,5)};
  Jacob(0,0)=1.0;
  Jacob(1,1)=1.0;
  Jacob(1,2)=-1.0;
  Jacob(2,2)=1.0;
  Jacob(2,3)=1.0;
  Jacob(3,4)=1.0;
  Jacob(4,3)=-K_eq;

  for(unsigned int i=0; i<Ca.rows(); ++i)
  {
   old_it<<Ca(i,step-1),H_piu(i,step-1),HCO3_meno(i,step-1),CO2(i,step-1),CaSiO3(i,step-1); //guess iniziale
   
   unsigned int max_iter=500;
   double tol=1.0e-14;
   double err=1; 
   Eigen::VectorXd dx(5);

      for(unsigned int iter=0; iter<max_iter and err>tol; ++iter)
       {
        compute_rhs(rhs, old_it, phi1(i), phi2(i), phi3(i), phi4(i), K_eq);///Calcolo F(x_k);
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
   }

}
  

void compute_rhs(Eigen::VectorXd& rhs,const Eigen::VectorXd old_it, double phi1, double phi2, double phi3, double phi4, double K_eq)
{

  double Ca(old_it(0)), H_piu(old_it(1)), HCO3_meno(old_it(2)), CO2(old_it(3)), CaSiO3(old_it(4));
  rhs<<Ca-phi1,H_piu-HCO3_meno-phi2,CO2+HCO3_meno-phi3,CaSiO3-phi4,H_piu*HCO3_meno-K_eq*CO2;

}


void compute_Jacob(Eigen::MatrixXd& J,const Eigen::VectorXd old_it)
{
     double HCO3_meno(old_it(2)),H_piu(old_it(1));
     J(4,1)=HCO3_meno;
     J(4,2)=H_piu;
}

