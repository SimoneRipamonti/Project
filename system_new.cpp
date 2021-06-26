//Loop:
//1) Calcolo dei vari rd.
//2) Funzione che mi risolve uno step temporale dei problemi di trasporto e reazione.
//3) Funzione che mi risolve il sistemino non lineare per trovare le concentrazioni effettive. 

#include <Eigen/Dense>
#include <cmath>

Eigen::MatrixXd::Zeros(Nx,Nt) Ca;
Eigen::MatrixXd::Zeros(Nx,Nt) H_piu;
Eigen::MatrixXd::Zeros(Nx,Nt) HCO3_meno;
Eigen::MatrixXd::Zeros(Nx,Nt) CO2;
Eigen::MatrixXd::Zeros(Nx,Nt) CaSiO3;
Eigen::MatrixXd::Zeros(Nx,Nt) SiO2;

Eigen::VectorXd::Zeros(Nx) phi1;
Eigen::VectorXd::Zeros(Nx) phi2;
Eigen::VectorXd::Zeros(Nx) phi3;
Eigen::VectorXd::Zeros(Nx) phi4;
Eigen::VectorXd::Zeros(Nx) phi5;


set_initial_cond(&Ca,&H_piu,&HCO3_meno,&CO2,&CaSiO3,&SiO2,Ca_0,H_piu_0,HCO3_meno_0,CO2_0,CaSiO3_0,SiO2_0);

//Setting for Transport and Reaction Equation
Matrix_F_piu F_p(Nx,Nx);
F_p.assemble_matrix(bc_cond,C_in,vel);

Matrix_F_meno F_m(Nx,Nx);
F_m.assemble_matrix(bc_cond,C_out,vel);

Matrix_C C(Nx,Nx);
C.assemble_matrix(phi,h);

const Eigen::MatrixXd M(1/dt*C.get_matrix()+F_p.get_matrix()-F_m.get_matrix());
const auto M_lu=M.fullPivLu();
  
const Eigen::MatrixXd rhs(Nx,Nx);
rhs=(1/dt*C.get_matrix());


const double const_r= area*rate_const*(std::exp(-E/(R*temperature)));

Eigen::VectorXd::Zeros(Nx) rd;

//Loop Temporale

for(unsigned int step=1; step<Nt; ++step)
{ 
  compute_phi(&phi1,&phi2,&phi3,&phi4,&phi5,Ca.col(step-1),H_piu.col(step-1),HCO3_meno.col(step-1),CO2.col(step-1),CaSiO3.col(step-1),SiO2.col(step-1)); //Calcolo le phi
  
  compute_rd(&rd,Ca.col(step-1),H_piu.col(step-1),SiO2.col(step-1),const_r,K_eq,n);//Calcolo i termini di reazione
  
  one_step_transport_reaction(&phi1,&phi2,&phi3,&phi4,&phi5,rd,M_lu,rhs); //Calcolo un passo della Reazione
  
  compute_concentration(&Ca,&H_piu,&HCO3,&CO2,&CaSiO3,&SiO2,phi1,phi2,phi3,phi4,phi5,step,K_eq); //Calcolo le Concentrazioni effettive
}

//Calcolo delle phi
compute_phi(Eigen::VectorXd& phi1,Eigen::VectorXd& phi2,Eigen::VectorXd& phi3,Eigen::VectorXd& phi4,Eigen::VectorXd& phi5,Eigen::VectorXd Ca,Eigen::VectorXd H_piu,Eigen::VectorXd HCO3_meno, Eigen::VectorXd CO2,Eigen::VectorXd CaSiO3,Eigen::VectorXd SiO2) 
{
phi1=Ca;
phi2=H_piu-HCO3_meno;
phi3=CO2+HCO3_meno;
phi4=CaSiO3;
phi5=SiO2;
}

//Calcolo Rd
compute_rd(&rd, const Eigen::VectorXd Ca, const Eigen::VectorXd H_piu, const Eigen::VectorXd SiO2,const double const_r, double K_eq, double n) //Calcolo il termine di reazione
{
const Eigen::VectorXd temp=const_r*H_piu.array().pow(n);
const Eigen::VectorXd omega=(Ca.cwiseProduct(SiO2)).cwiseQuotient(H_piu)/K_eq;
rd=temp.cwiseProduct(Eigen::VectorXd::Ones(Ca.size())-omega); 
}

//Un passo di reaction
one_step_transport_reaction(Eigen::VectorXd& phi1,Eigen::VectorXd& phi2,Eigen::VectorXd& phi3,Eigen::VectorXd& phi4,Eigen::VectorXd& phi5,const Eigen::VectorXd rd,const auto M_lu,const Eigen::MatrixXd rhs)
{
  transport_and_reaction(&phi1,M_lu,rhs,rd);
  transport_and_reaction(&phi2,M_lu,rhs,Eigen::Vector::Zeros(phi2.size());//no reaction
  transport_and_reaction(&phi3,M_lu,rhs,Eigen::Vector::Zeros(phi3.size());//no reaction
  transport_and_reaction(&phi4,M_lu,rhs,-rd);
  transport_and_reaction(&phi5,M_lu,rhs,rd);
}

//Trasporto e reazione
transport_and_reaction(Eigen::VectorXd& phi, const auto M_lu, const Eigen::MatrixXd &rhs, const Eigen::VectorXd rd)
{
  const Eigen::VectorXd temp=rhs*phi+rd;
  phi=M_lu.solve(temp);
}

//Risoluzione sistema non lineare
compute_concentration(Eigen::MatrixXd& Ca, Eigen::MatrixXd& H_piu, Eigen::MatrixXd& HCO3_meno, Eigen::MatrixXd& CO2, Eigen::MatrixXd& CaSiO3, Eigen::MatrixXd& SiO2, const Eigen::VectorXd phi1, const Eigen::VectorXd phi2, const Eigen::VectorXd phi3, const Eigen::VectorXd phi4, const Eigen::VectorXd phi5,unsigned int step,double K_eq)
{
//Initial guess::value at previous step

EigenVectorXd previous_iterate(6*Nx);
previous_iterate<<Ca.col(step-1),H_piu.col(step-1),HCO3_meno.col(step-1),CO2.col(step-1),CaSiO3.col(step-1),SiO2.col(step-1);


//Definizione rhs 
const Eigen::VectorXd K(Nx);
K.fill(K_eq);
Eigen::VectorXd rhs(phi1.size()+phi2.size()+phi3.size()+phi4.size()+phi5.size()+K.size());//rhs(6*Nx)
rhs<<phi1,phi2,phi3,phi4,phi5,K;

unsigne int max_iter=100;
double tol=1e-3;
double residuo=1;


//Definizione Jacobiana: è un super matricione?
Eigen::MatrixXd Jacob (6*Nx,6*Nx);
Eigen::MatrixXd sol (6*Nx,6*Nx);
for(unsigned int iter=0;iter<max_iter and res>tol;++iter)
{
  //Calcolo Jacobiana valutata nella guess
     
  Jacob<<Eigen::MatrixXd::Identity(Nx,Nx)...
  //sol=previous_iterate-Jacob^-1*f(x(k));
   
  //Calcolo valore concentrazioni della nuova iterata
  
  //Calcolo il residuo

}
  

}
