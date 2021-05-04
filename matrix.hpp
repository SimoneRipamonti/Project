#ifndef MATRIX_HH
#define MATRIX_HH
#include "parameters.hpp"
#include <string>
#include<Eigen/Dense>

class AbstractMatrix
{
public:

AbstractMatrix(const  std::string & filenamedata);

//AbstractMatrix(const AbstractMatrix &p);

//AbstractMatrix(const string &filename);

virtual void set_matrix()=0;

virtual void set_BC()=0;

virtual void set_rhs()=0;

virtual ~AbstractMatrix()=default;

protected:
unsigned int row;
unsigned int col;
Eigen::MatrixXd m;
Data data;
Eigen::VectorXd rhs;

};

class Matrix_A: public AbstractMatrix
{
public:
Matrix_A(const std::string &);
void set_matrix() override;
void set_BC() override; 
void set_rhs() override;
~Matrix_A(){};

};

class Matrix_B:public AbstractMatrix
{
public:
Matrix_B(const std::string &);
void set_matrix() override;
void set_BC() override;
void set_rhs() override;
~Matrix_B(){};
};

class Matrix_C:public AbstractMatrix
{
public:
  Matrix_C(const std::string &);
  void set_matrix() override;
  void set_BC() override;
  void set_rhs() override;
  ~Matrix_C(){};
};

class Matrix_F_piu:public AbstractMatrix
{
public:
   Matrix_F_piu(const std::string &);
   void set_matrix() override;
   void set_BC() override;
   void set_rhs() override;
   void set_velocity(const Eigen::VectorXd &vel);
   ~Matrix_F_piu(){};
   
   private:
   Eigen::VectorXd velocity;
   
};

class Matrix_F_meno:public AbstractMatrix
{
public:
   Matrix_F_meno(const std::string &);
   void set_matrix() override;
   void set_BC() override;
   void set_rhs() override;
   void set_velocity(const Eigen::VectorXd &vel);
   ~Matrix_F_meno(){};
   
   private:
   Eigen::VectorXd velocity;
};






class Matrix_
#endif

