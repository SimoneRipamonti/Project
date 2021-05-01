#ifndef MATRIX_HH
#define MATRIX_HH
#include "parameters.hpp"
#include <string>
#include<Eigen/Dense>

class AbstractMatrix
{
public:

AbstractMatrix(unsigned int r,unsigned int c,const  std::string & filenamedata);

AbstractMatrix(const AbstractMatrix &p);

//AbstractMatrix(const string &filename);

virtual void get_matrix()==0;

virtual void set_BC()==0;

virtual void set_rhs==0;

virtual ~AbstractMatrix()=default;

protected:
unsigned int rows;
unsigned int col;
Eigen::MatrixXd m;
Data data;
Eigen::VectorXd rhs;

};

class Matrix_A: public AbstractMatrix
{
public:
Matrix_A();
void get_matrix() override;
void set_BC() override; 
~Matrix_A(){};

};

class Matrix_B:public AbstractMatrix
{
public:
Matrix_B();
void get_matrix();
void set_BC() override;
~Matrix_B(){};

};

class Matrix_C:public AbstractMatrix
{
public:
void get_matrix();
~Matrix_C(){}
};
