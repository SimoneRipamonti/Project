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


Matrix_A::Matrix_A(unsigned int row,unsigned int col,const muparser_fun &per,double h_, double mu_,const std::string &inf,const std::string &out, double in_,double out_):AbstractMatrix(row,col),K(per),h(h_),mu(mu_),inflow(inf),outflow(out){

if(inflow=="Flow")
     Q_in=in_;
  else
     P_in=in_;
  if(outflow=="Flow")
    Q_out=out_;
  else
    P_out=out_;
}


void Matrix_A::set_matrix()
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
  else if(inflow=="Pressure")
     rhs(0)+=P_in;
  else
    throw std::invalid_argument("Invalid argument: wrong inflow boundary cond ");
      
  
  if(outflow=="Flow")
    {
      for(unsigned int i=0;i<col;++i)
           m(row-1,i)=0.;
      m(row-1,col-1)=1.;
      rhs(row-1)+=Q_out;
    }
  else if(outflow=="Pressure")
       rhs(row-1)-=P_out;
  else
    throw std::invalid_argument("Invalid argument: wrong outflow boundary cond");
      
}

void Matrix_A::set_rhs()
{}


Matrix_B::Matrix_B(unsigned int row, unsigned int col,const std::string inf_,const std::string out_,const muparser_fun &f,double h_):AbstractMatrix(row,col),inflow(inf_),outflow(out_),source(f),h(h_){}

void Matrix_B::set_matrix()
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


Matrix_C::Matrix_C(unsigned int row, unsigned int col,const std::string &bc,const muparser_fun &por_,double h_,double cond_):AbstractMatrix(row,col),bc_cond(bc),por(por_),h(h_)
{ 
 if(bc=="In")
   c_in=cond_;//è il flusso entrante c*vel
 else 
   c_out=cond_;//è il flusso uscente c*vel}
}

void Matrix_C::set_matrix()
{ 
 for(unsigned int i=0;i<row;++i)
      m(i,i)=por(0.5*h+i*h)*h;
}

void Matrix_C::set_BC()
{
  if(bc_cond=="In")
      {//m(0,0)=1.;
       rhs(0)+=c_in*por(0.5*h);
      }
  else if(bc_cond=="Out")
      {//m(row-1,row-1)=1.;
       rhs(row-1)-=c_out*(por((row-1)*h+0.5));
      }
  else 
    throw std::invalid_argument("Invalid argument: wrong boundary cond type");
}

void Matrix_C::set_rhs()
{}

Matrix_F_piu::Matrix_F_piu(unsigned int row, unsigned int col,const std::string &bc,const Eigen::VectorXd &vel):AbstractMatrix(row,col),bc_cond(bc),velocity(vel)
{}

void Matrix_F_piu::set_matrix()
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

void Matrix_F_piu::set_BC()
{
   if(bc_cond=="In")
       {//{m(0,0)=0.;
        //m(0,1)=0.;
       }
   else
        {//m(row-1,col-1)=0.;
 }
}

void Matrix_F_piu::set_rhs(){}

Matrix_F_meno::Matrix_F_meno(unsigned int row, unsigned int col,const std::string &bc,const Eigen::VectorXd &vel):AbstractMatrix(row,col),bc_cond(bc),velocity(vel)
{}

void Matrix_F_meno::set_matrix()
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

void Matrix_F_meno::set_BC(){
 
   if(bc_cond=="In")
       {//m(0,0)=0.;
       }
   else
       {
        //m(row-1,col-1)=0.; 
        //m(row-1,col-2)=0.;
       }
}

void Matrix_F_meno::set_rhs(){}



