/** 
 *  \file   matrix.hpp 
 *  \brief  discretized equations' matrices 
 *  \author Simone Ripamonti
 *  \date   2021-05 
 ***********************************************/

#ifndef MATRIX_HH
#define MATRIX_HH
#include "muparser_fun.hpp"
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>

/**
  *\brief Definition of an Abstract class for the matrices of the Darcy and Transport System.
  *
  * This class contain different virtual functions that will be overwritten by its derived classes.
  */
class AbstractMatrix
{
public:
     /*! \brief Constructor for the Abstract Class.
      *
      * It takes as input only the rows' and columns' number as  unsigned int.  
      */
    AbstractMatrix(unsigned int row, unsigned int col);

    /**
     * Function that defines the specific algebraic operator (the matrix) in question.
     */

    virtual void define_matrix()=0;

    virtual void set_BC()=0; /*!<Function that imposes the BC on the matrix*/

    virtual void set_rhs()=0;/*!<Function that sets the right-hand side corresponding to the problem*/
  
    Eigen::SparseMatrix<double>& get_matrix(); /*!<Function that gives in return the sparse matrix*/

    Eigen::VectorXd& get_rhs();/*!<Function that gives in return the right-hand side*/

    void print_m() const;/*!<Function that prints the matrix*/
    
    /**
     *Virtual destructor.
     */
    virtual ~AbstractMatrix()=default;

protected:
    unsigned int row; /*!<Number of rows*/
    unsigned int col;/*!<Number of columns*/
    Eigen::SparseMatrix<double> m;/*!<Sparse matrix that represent the algebraic operator*/
    Eigen::VectorXd rhs;/*!<Right-hand side of the algebraic operator*/


};



/*!
  *\brief Matrix A is the mass velocity matrix of the Darcy problem
*/
  
class Matrix_A: public AbstractMatrix
{
public:

    Matrix_A(unsigned int row,unsigned int col);
/*!
*Function that takes as input and sets the specific data necessary for assemblying matrix A
*\param per is the soil permeability passed as a muparser_fun since it changes in space
*\param h is the spatial step
*\param mu is the fluid viscosity
*\param inf is the type of in-section boundary condition data (velocity or pressure) 
*\param outf is the type of out-section boundary condition data (velocity or pressure)
*\param p_in is the in-pressure condition 
*\param p_out is the out-pressure condition
*\param q_in  is the inflow bc
*\param q_out is the outflow bc
*/
    void set_data(const muparser_fun &per,double h, double mu,const std::string &inf,const std::string &outf, double p_in, double p_out, double q_in, double q_out); 

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

/*!
 *Function that calls the above functions define_matrix(),set_BC(), set_rhs() in order to assemble all the matrix A
 */
    void assemble_matrix(const muparser_fun &per,double h, double mu,const std::string &inf,const std::string &outf, double p_in, double p_out, double q_in, double q_out);

    ~Matrix_A() {};

private:
    muparser_fun K; /*!<Permeablity of the soil passed as a muparser_fun since it is a function of the space*/
    double h; /*!<Spatial step*/
    double mu; /*!< Viscosity*/
    std::string inflow;/*!<Type of in-section boundary condition data (pressure or velocity)*/
    std::string outflow;/*!<Type of out-section data (pressure or velocity)*/
    double Q_in=0; /*!<Inflow data*/
    double Q_out=0;/*!<Outflow data*/
    double P_in=0;/*!<In-pressure data*/
    double P_out=0;/*!<Out-pressure data*/
};


/*!
  *\brief Matrix B is the saddle_point matrix of the Darcy_System.
*/
  
class Matrix_B:public AbstractMatrix
{
public:
    Matrix_B(unsigned int row, unsigned int col);
/*!
*Function that takes as input and sets the specific data necessary for assemblying matrix B
*\param inf is the type of in-section boundary condition data (velocity or pressure) 
*\param out is the type of out-section boundary condition data (velocity or pressure)
*\param f is the external source of the Darcy system's second equation
*\param h is the spatial step
*/
    void set_data(const std::string inf,const std::string out,const muparser_fun &f,double h);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const std::string inf,const std::string out,const muparser_fun &f,double h);/*!<Function that calls the above functions define_matrix(),set_BC(), set_rhs() in order to assemble all the matrix B*/

    ~Matrix_B() {};

private:
    std::string inflow;/*!<Type of in-section data (pressure or velocity)*/
    std::string outflow;/*!<Type of out-section data (pressure or velocity)*/
    muparser_fun source;/*!<External source f*/
    double h;/*!< Spatial step */
};

/*!
  *\brief Matrix C is the mass matrix for the Transport problem.
*/
  
class Matrix_C:public AbstractMatrix
{
public:
    Matrix_C(unsigned int row, unsigned int col);

   /*!
*Function that takes as input and sets the specific data necessary for assemblying matrix C
*\param por is the soil porosity
*\param h is the spatial step
*/   
   void set_data(const muparser_fun por,double h);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const muparser_fun por,double h);/*!<Function that calls the above functions define_matrix(),set_BC(), set_rhs() in order to assemble all the matrix C*/

    ~Matrix_C() {};

private:
    muparser_fun por;/*!<Soil porosity*/
    double h; /*!<Spatial step*/
};


/*!
  *\brief Matrix F_piu is the part of the Upwind Matrix which treats the left node of the cells.
*/
class Matrix_F_piu:public AbstractMatrix
{
public:
    Matrix_F_piu(unsigned int row, unsigned int col);

  /*!
*Function that takes as input and sets the specific data necessary for assemblying matrix F_piu
*\param bc type of boundary condition (in or out)
*\param c_bc tracer bc value
*\param vel transport fluid velocity
*/ 

    void set_data(const std::string &bc, double c_bc, const Eigen::VectorXd& vel);

    void define_matrix() override;

    void set_BC() override; 
   
    void set_rhs() override;

    void assemble_matrix(const std::string &bc, double c_bc, const Eigen::VectorXd& vel);/*!<Function that calls the above functions define_matrix(),set_BC(), set_rhs() in order to assemble all the matrix F_piu*/

    
    ~Matrix_F_piu() {};

private:
    std::string bc_cond;/*!<Type of bc condition (Inflow or Outflow)*/
    Eigen::VectorXd velocity;/*!<Transport fluid velocity*/
    double c_bc;/*!<tracer bc value*/
};

/*!
  *\brief Matrix F_meno is the part of the Upwind Matrix which treats the right node of the cells.
*/
class Matrix_F_meno:public AbstractMatrix
{
public:
    Matrix_F_meno(unsigned int row, unsigned col);

     /*!
*Function that takes as input and sets the specific data necessary for assemblying matrix F_meno
*\param bc type of boundary condition (in or out)
*\param c_bc tracer bc value
*\param vel transport fluid velocity
*/ 
    void set_data(const std::string &bc, double c_bc, const Eigen::VectorXd &vel);
 
    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;

    void assemble_matrix(const std::string &bc, double c_bc, const Eigen::VectorXd &vel);/*!<Function that calls the above functions define_matrix(),set_BC(), set_rhs() in order to assemble all the matrix F_meno*/

    ~Matrix_F_meno() {};

private:
    std::string bc_cond;/*!<Type of bc condition (Inflow or Outflow)*/
    Eigen::VectorXd velocity;/*!<Transport fluid velocity*/
    double c_bc;/*!<tracer bc value*/
};

/*!
  *\brief Matrix R is the Reaction matrix.
*/
class Matrix_R:public AbstractMatrix
{
public:
    Matrix_R(unsigned int row, unsigned int col);

     /*!
*Function that takes as input and sets the specific data necessary for assemblying matrix R
*\param area is the reaction Area
*\param rate_const is the constant rate
*\param temperature is the temperature reaction
*\param R is the gas constant
*\param E is the activation energy
*\param ph is the water ph
*\param const_eq is the equilibrium reaction constant
*/ 

    void set_data(double area, double rate_const, double temperature, double R, double E, double ph,double const_eq, double h);

    void define_matrix() override;

    void set_BC() override;

    void set_rhs() override;
   
    void assemble_matrix(double area, double rate_const, double temperature, double R, double E, double ph,double const_eq, double h);/*!<Function that calls the above functions define_matrix(),set_BC(), set_rhs() in order to assemble all the matrix */

    void update(const Eigen::VectorXd &past_sol);/*!<Function that updates the Reaction matrix since we treat this part explicitely*/

    ~Matrix_R() {};

private:
    double react_const;/*!<reaction constant computed with data*/
    double ph;/*!<water ph*/
    double const_eq;/*!<Equilibrium reaction constant*/
    double h;/*!<Spatial step*/
};




#endif

