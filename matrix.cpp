#include "matrix.hpp"
#include <exception> 
#include <cmath>
AbstractMatrix::AbstractMatrix(const std::string &filenamedata):data(filenamedata){}

//AbstractMatrix::AbstractMatrix(const AbstractMatrix &p):rows(p.rows),col(p.col),m(p.rows,p.col),rhs(p.rhs),data(p.data){}

Matrix_A::Matrix_A(const std::string &filenamedata):AbstractMatrix(filenamedata){
rows=data.Nx+1;
col=data.Nx+1;
auto a=data.Nx+1;
m(a,a);
rhs(a);
}

Matrix_B::Matrix_B(const std::string &filenamedata):AbstractMatrix(filenamedata){
rows=data.Nx;
col=data.Nx;
auto a=data.Nx+1;
m(data.Nx,a);
rhs(data.Nx);
}

void Matrix_A::get_matrix()
{
 double h=static_cast<double>(data.domain_length)/(data.Nx);
 unsigned int Nx=data.Nx;
 auto K=data.K;
 double mu=data.mu;
 m(0,0)=1./3.*h*std::pow(K(h),-1)/mu;
 m(1,0)=1./6.*h*std::pow(K(h),-1)/mu;
 for(size_t i=1;i<rows-2;++i)
 {
   m(i,i-1)=1./6.*h*std::pow(K(i*h),-1)/mu;
   m(i,i)=1./3.*h*std::pow(K(i*h),-1)/mu+1./3.*h*std::pow(K((i+1)*h),-1)/mu;
   m(i,i+1)=1./6.*h*std::pow(K((i+1)*h),-1)/mu;
 }
 m(Nx-1,Nx-2)=1./6.*h*std::pow(K((Nx-2)*h),-1)/mu;
 m(Nx-1,Nx-1)=1./3.*h*std::pow(K((Nx-2)*h),-1)/mu;
}

void Matrix_A::set_BC()
{  
  if(data.BC_in=="Flow")
     {
      for(int i=0;i<m.cols();++i)
            m(0,i)=0.;
      m(0,0)=1.;
      rhs(0)+=data.Q_in;
     }
  else if(data.BC_in=="Pressure")
     {
      rhs(0)+=data.p_in;
     }
  else{
    throw std::invalid_argument("Invalid argument: wrong inflow boundary cond ");
      }
  
  if(data.BC_out=="Flow")
    {
      for(int i=0;i<m.cols();++i)
           m(rows-1,i)=0.;
      m(rows-1,rows-1)=1.;
      rhs(rows-1)+=data.Q_out;
    }
  else if(data.BC_out=="Pressure")
    {
     rhs(rows-1)-=data.p_out;
    }
  else{
    throw std::invalid_argument("Invalid argument: wrong outflow boundary cond");
      }
}

void Matrix_A::set_rhs()
{
   for(int i=0;i<rhs.size();++i)
         rhs(i)=0.;
}

void Matrix_B::get_matrix()
{
  for(unsigned int i=0;i<rows-1;++i)
    {
      m(i,i)=-1;
      m(i,i+1)=1;
    }
  m(rows-1,col-1)=-1;
}


void Matrix_B::set_BC()
{ 
  if(data.BC_in=="Flow")
      {for(int i=0;i<m.cols();++i)
            m(0,i)=0.;}
  if(data.BC_out=="Flow")
     {for(int i=0;i<m.cols();++i)
           m(rows-1,i)=0.;}
}

void Matrix_B::set_rhs()
{
  auto source=data.f;
  double h=static_cast<double>(data.domain_length)/(data.Nx);
  for(size_t i=0;i<rows-2;++i)
    rhs(i)=source(i*h);
}
