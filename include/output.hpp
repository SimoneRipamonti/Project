#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"


/** \addtogroup Output_Function
 *\brief Functions that print on a csv file the solutions of spatial-temporal problems (the transport and reaction ones)
 *  @{
*/
/*!
 *Function that saves the solution in a .csv file the following way: each row is a spatial position, each column is a time instant
*\param title is the file title 
*\param value1 is a reference to the matrix where we have stored our solution
*\param L is the domain length
*\param Nx is the number of spatial steps
*\param T is the temporal domain
*\param Nt is the number of temporal steps
*/
void output_results_fixed_time(const std::string& title, const Matrix_full &value1, double L, unsigned int Nx, double T,unsigned int Nt);


/*!
 *Function that saves the solution in a .csv file the following way: each row is a time instant, each column is a spatial position 
*\param title is the file title 
*\param value1 is a reference to the matrix where we have stored our solution
*\param L is the domain length
*\param Nx is the number of spatial steps
*\param T is the temporal domain
*\param Nt is the number of temporal steps
*/
void output_results_fixed_space(const std::string& title, const Matrix_full &value1, double L, unsigned int Nx, double T,unsigned int Nt);

void output_example(const Vector &sol,unsigned int Nx, double L);
/** @}*/



#endif

