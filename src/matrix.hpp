#ifndef MATRIX_HH
#define MATRIX_HH
#include "muparser_fun.hpp"
#include <string>
#include <Eigen/Dense>

//Definition of an Abstract class for the matrices of the Darcy and Transport System
class AbstractMatrix
{
public:

    AbstractMatrix(unsigned int row, unsigned int col);

    virtual void define_matrix()=0;

    virtual void set_BC()=0;

    virtual void set_rhs()=0;

    Eigen::MatrixXd& get_matrix();

    Eigen::VectorXd& get_rhs();

    void print_m() const;
    virtual ~AbstractMatrix()=default;

protected:
    unsigned int row;
    unsigned int col;
    Eigen::MatrixXd m;
    Eigen::VectorXd rhs;

};


//Matrix_A is the mass matrix of the Darcy Problem
class Matrix_A: public AbstractMatrix
{
public:

    Matrix_A(unsigned int row,unsigned int col);

//set_data function takes as input the spcific data that we need for assembling the matrix in question
    void set_data(const muparser_fun &per,double h, double mu,const std::string &inf,const std::string &outf, double in,double out);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

//assemble_matrix simply call the above functions define_matrix(),set_BC(),set_rhs() in order to assemble the matrix considering the BC and the rhs
    void assemble_matrix(const muparser_fun &per,double h, double mu,const std::string &inf,const std::string &outf, double in,double out);

    ~Matrix_A() {};

private:
    muparser_fun K;
    double h;
    double mu;
    std::string inflow;
    std::string outflow;
    double Q_in=0;
    double Q_out=0;
    double P_in=0;
    double P_out=0;
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
    std::string inflow;
    std::string outflow;
    muparser_fun source;
    double h;
};


//Matrix_C is the mass matrix of the Transport_System
class Matrix_C:public AbstractMatrix
{
public:
    Matrix_C(unsigned int row, unsigned int col);

    void set_data(const std::string &bc,const muparser_fun por,double h,double cond);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const std::string &bc,const muparser_fun por,double h,double cond);

    ~Matrix_C() {};

private:
    std::string bc_cond;
    muparser_fun por;
    double h;
    double c_in=0;
    double c_out=0;
};


//Matrix F_piu it's the part of the Upwind Matrix which treats the left node of the cells
class Matrix_F_piu:public AbstractMatrix
{
public:
    Matrix_F_piu(unsigned int row, unsigned int col);

    void set_data(const std::string &bc,const Eigen::VectorXd &vel);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const std::string &bc,const Eigen::VectorXd &vel);

    ~Matrix_F_piu() {};

private:
    std::string bc_cond;
    Eigen::VectorXd velocity;
};

//Matrix F_meno it's the part of the Upwind Matrix which treats the right node of the cells
class Matrix_F_meno:public AbstractMatrix
{
public:
    Matrix_F_meno(unsigned int row, unsigned col);

    void set_data(const std::string &bc,const Eigen::VectorXd &vel);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const std::string &bc,const Eigen::VectorXd &vel);

    ~Matrix_F_meno() {};

private:
    std::string bc_cond;
    Eigen::VectorXd velocity;
};

class Matrix_R:public AbstractMatrix
{
public:
    Matrix_R(unsigned int row, unsigned int col);

    void set_data(double area, double rate_const, double temperature, double R, double E, double ph,double const_eq);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;
   
    void assemble_matrix(double area, double rate_const, double temperature, double R, double E, double ph,double const_eq);

    void update(const Eigen::VectorXd &past_sol);

    ~Matrix_R() {};

private:
    double react_const;
    double ph;
    double const_eq;
};




#endif

