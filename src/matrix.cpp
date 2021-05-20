#include "matrix.hpp"
#include <exception> 
#include <cmath>
AbstractMatrix::AbstractMatrix(unsigned int row_,unsigned int col_):row(row_),col(col_){
m=Eigen::MatrixXd::Zero(row_,col_);
rhs=Eigen::VectorXd::Zero(col_); 
}

Eigen::MatrixXd& AbstractMatrix::get_matrix(){
     return m;
}

Eigen::VectorXd& AbstractMatrix::get_rhs(){
     return rhs;
}

void AbstractMatrix::print_m(){
std::cout<<m<<std::endl;
}


Matrix_A::Matrix_A(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_){}

void Matrix_A::set_data(const muparser_fun &per_,double h_, double mu_,const std::string &inf_,const std::string &outf_, double in_,double out_)
{
 K=per_;
 h=h_;
 mu=mu_;
 inflow=inf_;
 outflow=outf_;
 
 if(inflow=="Flow")
     Q_in=in_;
 else if(inflow=="Pressure")
         P_in=in_;
      else
         throw std::invalid_argument("Invalid argument: wrong inflow boundary cond ");
 
 if(outflow=="Flow")
    Q_out=out_;
 else if(outflow=="Pressure")
        P_out=out_;
      else
        throw std::invalid_argument("Invalid argument: wrong inflow boundary cond ");
}

void Matrix_A::define_matrix()
{
 m(0,0)=1./3.*h*std::pow(K(h)/mu,-1);
 m(0,1)=1./6.*h*std::pow(K(h)/mu,-1);
 for(unsigned int i=1;i<row-1;++i)
 {
   m(i,i-1)=1./6.*h*std::pow(K(i*h)/mu,-1);
   m(i,i)=1./3.*h*std::pow(K(i*h)/mu,-1)+1./3.*h*std::pow(K((i+1)*h)/mu,-1);
   m(i,i+1)=1./6.*h*std::pow(K((i+1)*h)/mu,-1);
 }
 m(row-1,col-2)=1./6.*h*std::pow(K((row-1)*h)/mu,-1);
 m(row-1,col-1)=1./3.*h*std::pow(K((row-1)*h)/mu,-1);
}

void Matrix_A::set_BC()
{  
  if(inflow=="Flow")
     {
      for(unsigned int i=0;i<col;++i)
            m(0,i)=0.;
      m(0,0)=1.;
      rhs(0)+=Q_in;
     }
  else 
     rhs(0)+=P_in;
  
  if(outflow=="Flow")
    {
      for(unsigned int i=0;i<col;++i)
           m(row-1,i)=0.;
      m(row-1,col-1)=1.;
      rhs(row-1)+=Q_out;
    }
  else
    rhs(row-1)-=P_out;
}

void Matrix_A::set_rhs()
{}

void Matrix_A::assemble_matrix(const muparser_fun &per_,double h_, double mu_,const std::string &inf_,const std::string &outf_, double in_,double out_){
set_data(per_,h_,mu_,inf_,outf_,in_,out_);
define_matrix();
set_BC();
set_rhs();
}


Matrix_B::Matrix_B(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_){}

void Matrix_B::set_data(const std::string inf_,const std::string out_,const muparser_fun &f_,double h_)
{
inflow=inf_;
outflow=out_;
source=f_;
h=h_;
}

void Matrix_B::define_matrix()
{
  m(0,0)=1;
  for(unsigned int i=1;i<col;++i)
    {
      m(i,i)=1;
      m(i,i-1)=-1;
    }
  m(row-1,col-1)=-1;
}


void Matrix_B::set_BC()
{ 
if(inflow=="Flow")
     {for(unsigned int i=0;i<col;++i)
            m(0,i)=0.;}
if(outflow=="Flow")
     {for(unsigned int i=0;i<col;++i)
           m(row-1,i)=0.;}
}

void Matrix_B::set_rhs()
{
  for(float i=0.5;i<col;++i)
    rhs(i)=source(i*h)*h;
}

void Matrix_B::assemble_matrix(const std::string inf_,const std::string out_,const muparser_fun &f_,double h_){
set_data(inf_,out_,f_,h_);
define_matrix();
set_BC();
set_rhs();
}

Matrix_C::Matrix_C(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_){}

void Matrix_C::set_data(const std::string &bc_,const muparser_fun &por_,double h_,double cond_){
 bc_cond=bc_;
 por=por_;
 h=h_;
  if(bc_cond=="In")
    c_in=cond_;
  else if(bc_cond=="Out")
          c_out=cond_;
       else 
          throw std::invalid_argument("Invalid argument: wrong boundary cond type");
}

void Matrix_C::define_matrix()
{ 
 for(unsigned int i=0;i<row;++i)
      m(i,i)=por(0.5*h+i*h)*h;
}

void Matrix_C::set_BC()
{
  if(bc_cond=="In")
     rhs(0)+=c_in*por(0.5*h);
  else 
     rhs(row-1)-=c_out*(por((row-1)*h+0.5));
      
}

void Matrix_C::set_rhs()
{}

void Matrix_C::assemble_matrix(const std::string &bc_,const muparser_fun &por_,double h_,double cond_){
set_data(bc_,por_,h_,cond_);
define_matrix();
set_BC();
set_rhs();
}

Matrix_F_piu::Matrix_F_piu(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_){}

void Matrix_F_piu::set_data(const std::string &bc_,const Eigen::VectorXd &vel_)
{
bc_cond=bc_;
velocity=vel_;
}

void Matrix_F_piu::define_matrix()
{
  for(unsigned int i=0;i<row-1;++i)
      { 
         if(velocity(i+1)>0)
             m(i,i)=velocity(i+1);
          else
             m(i,i+1)=velocity(i+1);   
      }
      m(row-1,col-1)=velocity(row);
}

void Matrix_F_piu::set_BC(){}

void Matrix_F_piu::set_rhs(){}

void Matrix_F_piu::assemble_matrix(const std::string &bc_,const Eigen::VectorXd &vel_){
set_data(bc_,vel_);
define_matrix();
set_BC();
set_rhs();
}


Matrix_F_meno::Matrix_F_meno(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_){}

void Matrix_F_meno::set_data(const std::string &bc_,const Eigen::VectorXd &vel_)
{
bc_cond=bc_;
velocity=vel_;
}

void Matrix_F_meno::define_matrix()
{
 //m(0,0)=velocity(0);
 for(unsigned int i=1;i<row;++i)
    {
      if(velocity(i)>0)
          m(i,i-1)=velocity(i);
      else
          m(i,i)=velocity(i);
    }
}

void Matrix_F_meno::assemble_matrix(const std::string &bc_,const Eigen::VectorXd &vel_){
set_data(bc_,vel_);
define_matrix();
set_BC();
set_rhs();
}

void Matrix_F_meno::set_BC(){}

void Matrix_F_meno::set_rhs(){}



