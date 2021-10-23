#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"


/** \addtogroup Darcy_Output
 *\brief Functions that print on a csv file the solutions of Darcy problem 
 *  @{
*/
/*!
 *Function that saves the solution (velocity+pressure) of the Darcy system on two different csv files (one for pressure and one for velocity)
 *\param value is a refence to the vector which contains the values that we want to store
 *\param Nx is the number of spatial inervals
 *\param L is the domain lenght 
*/
void Darcy_output_results(Vector &value,unsigned int Nx,double L);

/*!
 *Function that saves the exact pressure values at the nodes point in a csv file
 *\param value1 is a refence to the vector which contains the values that we want to store
 *\param Nx is the number of spatial inervals
 *\param L is the domain lenght 
*/
void pressure_exact_result(Vector &value1, unsigned int Nx, double L);

/*!
 *Function that saves the exact velocity values at the nodes point in a csv file
 *\param value1 is a refence to the vector which contains the values that we want to store
 *\param Nx is the number of spatial inervals
 *\param L is the domain lenght 
*/
void velocity_exact_result(Vector &value1, unsigned int Nx, double L);

/*!
 *Function that saves the error values in a csv file
 *\param value1 is a refence to the vector which contains the values that we want to store
 *\param Nx is a reference to the vector which contains all the spatial steps used for the convergence analyses  
*/
void output_error(Vector &value1, Eigen::VectorXi &Nx, const std::string &name);

/** @}*/


#endif

