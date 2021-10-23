#include "matrix.hpp"
#include <exception>
#include <cmath>
#include <algorithm>

///AbstractMatrix constructor
AbstractMatrix::AbstractMatrix(unsigned int row_,unsigned int col_):row(row_),col(col_)
{
    m.resize(row_,col_);
    rhs=Eigen::VectorXd::Zero(col_);
}

/**
 *Getters
 */
Eigen::SparseMatrix<double>& AbstractMatrix::get_matrix() 
{
    return m;
}

Eigen::VectorXd& AbstractMatrix::get_rhs() 
{
    return rhs;
}

//Print the matrix
void AbstractMatrix::print_m() const
{
    std::cout<<m<<std::endl;
}

//Matrix_A constructor
Matrix_A::Matrix_A(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_) {}

//set_data gives to the matrix object the physical and geometrical parameters which it needs to assemble the matrix
void Matrix_A::set_data(const muparser_fun &per_,double h_, double mu_,const std::string &inf_,const std::string &outf_, double p_in_, double p_out_, double q_in_, double q_out_)
{
    K=per_;
    h=h_;
    mu=mu_;
    inflow=inf_;
    outflow=outf_;

//We have two possible inflow/outflow condition: or the one related to the velocity or the one related to the pressure
    if(inflow=="Flow")
        Q_in=q_in_;
    else if(inflow=="Pressure")
        P_in=p_in_;
    else
        throw std::invalid_argument("Invalid argument: wrong inflow boundary cond ");

    if(outflow=="Flow")
        Q_out=q_out_;
    else if(outflow=="Pressure")
        P_out=p_out_;
    else
        throw std::invalid_argument("Invalid argument: wrong inflow boundary cond ");
}

//Definition of the mass matrix A considering P1 Finite element for the Velocity and P0 for the pressure (They are stable for the 1D case)
void Matrix_A::define_matrix()
{
    m.insert(0,0)=1./3.*h*std::pow(K(h)/mu,-1);
    m.insert(0,1)=1./6.*h*std::pow(K(h)/mu,-1);
    for(unsigned int i=1; i<row-1; ++i)
    {
        m.insert(i,i-1)=1./6.*h*std::pow(K(i*h)/mu,-1);
        m.insert(i,i)=1./3.*h*std::pow(K(i*h)/mu,-1)+1./3.*h*std::pow(K((i+1)*h)/mu,-1);
        m.insert(i,i+1)=1./6.*h*std::pow(K((i+1)*h)/mu,-1);
    }
    m.coeffRef(row-1,col-2)=1./6.*h*std::pow(K((row-1)*h)/mu,-1);
    m.coeffRef(row-1,col-1)=1./3.*h*std::pow(K((row-1)*h)/mu,-1);
}

//set_BC() function change the first/last row of the matrix A considering what kind of BC are set
void Matrix_A::set_BC()
{
    if(inflow=="Flow")//if we have Dirichlet condition for the velocity we have to impose it in a strong way
    {
        m.coeffRef(0,1)=0.;
        m.coeffRef(0,0)=1.;
        rhs(0)+=Q_in;
    }
    else //if we have Dirichlet condition for the pressure we impose it in a weak way since we have the weak term in the equation
        rhs(0)+=P_in;

    if(outflow=="Flow")
    {

        m.coeffRef(row-1,col-2)=0.;
        m.coeffRef(row-1,col-1)=1.;
        rhs(row-1)+=Q_out;
    }
    else
        rhs(row-1)-=P_out;
}

//There is not any source term in the first equation of the Darcy system for the assumptions made
void Matrix_A::set_rhs()
{}


void Matrix_A::assemble_matrix(const muparser_fun &per_,double h_, double mu_,const std::string &inf_,const std::string &outf_, double p_in_, double p_out_, double q_in_, double q_out_)
{
    set_data(per_,h_,mu_,inf_,outf_,p_in_,p_out_,q_in_,q_out_);
    define_matrix();
    set_BC();
    set_rhs();
}

//Matrix_B constructor
Matrix_B::Matrix_B(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_) {}

//set_data gives to the matrix object the physical and geometrical parameters which it needs to assemble the matrix
void Matrix_B::set_data(const std::string inf_,const std::string out_,const muparser_fun &f_,double h_)
{
    inflow=inf_;
    outflow=out_;
    source=f_;
    h=h_;
}

//Definition of the B matrix with P1-P0 element
void Matrix_B::define_matrix()
{
    m.insert(0,0)=1;
    for(unsigned int i=1; i<col; ++i)
    {
        m.insert(i,i)=1;
        m.insert(i,i-1)=-1;
    }
    m.insert(row-1,col-1)=-1;
}

//set_BC() function change the first/last row of the matrix B considering what kind of BC we have
void Matrix_B::set_BC()
{
    if(inflow=="Flow")
        m.coeffRef(0,0)=0.;

    if(outflow=="Flow")
        m.coeffRef(row-1,col-1)=0.;

//If we have in Inflow/outflow pressure condition there is nothing to be made
}

//Here the source term that carachterize the second equation of the Darcy System is introduced (the one that is equal to the divergence of the velocity)
void Matrix_B::set_rhs()
{
    for(float i=0.5; i<col; ++i)
        rhs(i)=source(i*h)*h;
}

void Matrix_B::assemble_matrix(const std::string inf_,const std::string out_,const muparser_fun &f_,double h_)
{
    set_data(inf_,out_,f_,h_);
    define_matrix();
    set_BC();
    set_rhs();
}

//Matrix_C constructor
Matrix_C::Matrix_C(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_) {}

//set_data gives to the matrix object the physical and geometrical parameters which it needs to assemble the matrix

void Matrix_C::set_data(const muparser_fun phi_,double h_)
{
    phi=phi_;
    h=h_;
}


//Definition of Matrix_C considering that we use P1 element
void Matrix_C::define_matrix()
{
    for(unsigned int i=0; i<row; ++i)
        m.insert(i,i)=phi(0.5*h+i*h)*h;//we include in the mass matrix of the transport equation (Matrix_C) already the product with the spatial step h and the porosity

}

//The BC of the Transport problem are treated in the F matrices
void Matrix_C::set_BC()
{}

//There is not no source term of the Transport Equation
void Matrix_C::set_rhs()
{}

void Matrix_C::assemble_matrix(const muparser_fun phi_,double h_)
{
    set_data(phi_,h_);
    define_matrix();
    set_BC();
    set_rhs();
}

//Matrix_F_piu constructor
Matrix_F_piu::Matrix_F_piu(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_) {}

//set_data gives to the matrix object the physical and geometrical parameters which it needs to assemble the matrix
void Matrix_F_piu::set_data(const std::string &bc_, double c_bc_, const Vector& vel_)
{
    bc_cond=bc_;
    velocity=vel_;
    c_bc=c_bc_;
}

//Defintion of matrix F_piu as the left part of the upwind scheme with P1 element
void Matrix_F_piu::define_matrix()
{
    for(unsigned int i=0; i<row-1; ++i)
    {
        if(velocity(i+1)>0)
            m.insert(i,i)=velocity(i+1);
        else
            m.insert(i,i+1)=velocity(i+1);
    }

    if(velocity(row)>0)
        m.insert(row-1,col-1)=velocity(row);
}

void Matrix_F_piu::set_BC() {}

//Here we treat the BC of the transport problem working on the rhs of the matrix
void Matrix_F_piu::set_rhs()
{

    if(bc_cond=="In")
        rhs(0)+=velocity(1)*c_bc;

}

void Matrix_F_piu::assemble_matrix(const std::string &bc_, double c_bc_, const Vector& vel_)
{
    set_data(bc_,c_bc_,vel_);
    define_matrix();
    set_BC();
    set_rhs();
}


//Matrix_F_meno constructor
Matrix_F_meno::Matrix_F_meno(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_) {}

//set_data gives to the matrix object the physical and geometrical parameters which it needs to assemble the matrix
void Matrix_F_meno::set_data(const std::string &bc_, double c_bc_, const Vector &vel_)
{
    bc_cond=bc_;
    velocity=vel_;
    c_bc=c_bc_;
}

//Defintion of matrix F_meno as the right part of the upwind scheme with P1 element
void Matrix_F_meno::define_matrix()
{
    if(velocity(0)<0)
        m.insert(0,0)=velocity(0);
    for(unsigned int i=1; i<row; ++i)
    {
        if(velocity(i)>0)
            m.insert(i,i-1)=velocity(i);
        else
            m.insert(i,i)=velocity(i);
    }
}

void Matrix_F_meno::assemble_matrix(const std::string &bc_,double c_bc_, const Vector &vel_)
{
    set_data(bc_,c_bc_,vel_);
    define_matrix();
    set_BC();
    set_rhs();
}

void Matrix_F_meno::set_BC() {}

//Here we treat the BC of the transport problem working on the rhs of the matrix
void Matrix_F_meno::set_rhs()
{

    if(bc_cond=="Out")
        rhs(row-1)+=velocity(row-1)*c_bc;

}



//Matrix_R
Matrix_R::Matrix_R(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_) {}

//set_data gives to the matrix object the physical and geometrical parameters which it needs to assemble the matrix
void Matrix_R::set_data(double area, double rate_const, double temperature, double R, double E, double ph_, double const_eq_, double h_, const muparser_fun phi_)
{
    //react_const=area*rate_const*(std::exp(-E/(R*temperature)))*std::pow(10,-ph);
    react_const=area*rate_const*(std::exp(-E/(R*temperature)));
    ph=ph_;
    const_eq=const_eq_;
    h=h_;
    phi=phi_;
}

//Defintion of matrix F_piu as the right part of the upwind scheme with P1 elemen
void Matrix_R::define_matrix() {}

void Matrix_R::set_BC() {}

void Matrix_R::set_rhs() {}

void Matrix_R::assemble_matrix(double area, double rate_const, double temperature, double R, double E, double ph,double const_eq, double h, const muparser_fun phi)
{
    set_data(area,rate_const,temperature,R,E,ph,const_eq,h,phi);
    define_matrix();
    set_BC();
    set_rhs();
}


//Here the reaction matrix is update with the past_solution having in mind a particular square reaction rate
void Matrix_R::update(const Eigen::VectorXd &past_sol)
{

    const Eigen::VectorXd p=past_sol.array().pow(2)/(const_eq*std::pow(10,-2*ph));
    for(unsigned int i=0; i<rhs.size(); ++i)
        rhs(i)=h*phi(h/2+i*h)*std::max(react_const*(1.0-p(i)),0.0);

}

//Matrix_A constructor
A_example::A_example(unsigned int row_,unsigned int col_):AbstractMatrix(row_,col_) {}

//set_data gives to the matrix object the physical and geometrical parameters which it needs to assemble the matrix
void A_example::set_data(const muparser_fun &source_, unsigned Nx, double L)
{
   source=source_;
   h=L/Nx;
}

//Definition of the mass matrix A considering P1 Finite element for the Velocity and P0 for the pressure (They are stable for the 1D case)
void A_example::define_matrix()
{
    m.insert(0,0)=1.;
    m.insert(0,1)=0.;
    for(unsigned int i=1; i<row-1; ++i)
    {
        m.insert(i,i-1)=-1./std::pow(h,2);
        m.insert(i,i)=2./std::pow(h,2);
        m.insert(i,i+1)=-1./std::pow(h,2);
    }
    m.coeffRef(row-1,col-2)=0.;
    m.coeffRef(row-1,col-1)=1.;
}

//set_BC() function change the first/last row of the matrix A considering what kind of BC are set
void A_example::set_BC()
{   
    rhs(0)=0.;
    rhs(row-1)=0.;
}

//There is not any source term in the first equation of the Darcy system for the assumptions made
void A_example::set_rhs()
{{
    for(unsigned i=0; i<col; ++i)
        rhs(i)=source(i*h);
}
}


void A_example::assemble_matrix(const muparser_fun &source_, unsigned Nx, double L)
{
    set_data(source_, Nx, L);
    define_matrix();
    set_BC();
    set_rhs();
}







