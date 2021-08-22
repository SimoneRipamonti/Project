#ifndef MATRIX_HH
#define MATRIX_HH
#include "muparser_fun.hpp"
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>

//Definition of an Abstract class for the matrices of the Darcy and Transport System
class AbstractMatrix
{
public:

    AbstractMatrix(unsigned int row, unsigned int col);

    virtual void define_matrix()=0;

    virtual void set_BC()=0;

    virtual void set_rhs()=0;

    Eigen::SparseMatrix<double>& get_matrix();

    Eigen::VectorXd& get_rhs();

    void print_m() const;
    virtual ~AbstractMatrix()=default;

protected:
    unsigned int row;
    unsigned int col;
    Eigen::SparseMatrix<double> m;
    Eigen::VectorXd rhs;

};


//Matrix_A is the mass matrix of the Darcy Problem
class Matrix_A: public AbstractMatrix
{
public:

    Matrix_A(unsigned int row,unsigned int col);

//set_data function takes as input the spcific data that we need for assembling the matrix in question
    void set_data(const muparser_fun &per,double h, double mu,const std::string &inf,const std::string &outf, double p_in, double p_out, double q_in, double q_out);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

//assemble_matrix simply call the above functions define_matrix(),set_BC(),set_rhs() in order to assemble the matrix considering the BC and the rhs
    void assemble_matrix(const muparser_fun &per,double h, double mu,const std::string &inf,const std::string &outf, double p_in, double p_out, double q_in, double q_out);

    ~Matrix_A() {};

private:
    muparser_fun K; //permeablity of the soil
    double h; //spatial step
    double mu; //viscosity
    std::string inflow; //Type of in-section data (pressure or velocity)
    std::string outflow; ////Type of out-section data (pressure or velocity)
    double Q_in=0; //Inflow data
    double Q_out=0;//Outflow data
    double P_in=0;//In pressure data
    double P_out=0;//Out pressure data
};

//Matrix_B is the saddle_point matrix of the Darcy_System
class Matrix_B:public AbstractMatrix
{
public:
    Matrix_B(unsigned int row, unsigned int col);

    void set_data(const std::string inf,const std::string out,const muparser_fun &f,double h);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const std::string inf,const std::string out,const muparser_fun &f,double h);

    ~Matrix_B() {};

private:
    std::string inflow;//Type of in-section data (pressure or velocity)
    std::string outflow;//Type of out-section data (pressure or velocity)
    muparser_fun source;//External source
    double h;//Spatial step
};


//Matrix_C is the mass matrix of the Transport_System
class Matrix_C:public AbstractMatrix
{
public:
    Matrix_C(unsigned int row, unsigned int col);

    void set_data(const muparser_fun por,double h);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const muparser_fun por,double h);

    ~Matrix_C() {};

private:
    muparser_fun por;//porosity of the soil
    double h; //spatial step
};


//Matrix F_piu it's the part of the Upwind Matrix which treats the left node of the cells
class Matrix_F_piu:public AbstractMatrix
{
public:
    Matrix_F_piu(unsigned int row, unsigned int col);

    void set_data(const std::string &bc, double c_bc, const Eigen::VectorXd& vel);

    void define_matrix() override;

    void set_BC() override; 
   
    void set_rhs() override;

    void assemble_matrix(const std::string &bc, double c_bc, const Eigen::VectorXd& vel);
    
    ~Matrix_F_piu() {};

private:
    std::string bc_cond;//Type of bc condition (Inflow or Outflow)
    Eigen::VectorXd velocity;//Velocity
    double c_bc;
};

//Matrix F_meno it's the part of the Upwind Matrix which treats the right node of the cells
class Matrix_F_meno:public AbstractMatrix
{
public:
    Matrix_F_meno(unsigned int row, unsigned col);

    void set_data(const std::string &bc, double c_bc, const Eigen::VectorXd &vel);
 
    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const std::string &bc, double c_bc, const Eigen::VectorXd &vel);

    ~Matrix_F_meno() {};

private:
    std::string bc_cond;
    Eigen::VectorXd velocity;
    double c_bc;
};


//Reaction Matrix
class Matrix_R:public AbstractMatrix
{
public:
    Matrix_R(unsigned int row, unsigned int col);

    void set_data(double area, double rate_const, double temperature, double R, double E, double ph,double const_eq, double h);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;
   
    void assemble_matrix(double area, double rate_const, double temperature, double R, double E, double ph,double const_eq, double h);

    void update(const Eigen::VectorXd &past_sol);//This fuction updare the Reaction matrix since we treat this part explicitely

    ~Matrix_R() {};

private:
    double react_const;//constant rate
    double ph;//ph of the undergound water
    double const_eq;//eq constant of the reaction
    double h;//spatial step
};




#endif

