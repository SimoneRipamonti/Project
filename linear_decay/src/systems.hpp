#ifndef SYSTEMS_HH
#define SYSTEMS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"



/*!
*Definition of the Transport system solved with an explicit upwind scheme
*/
void Transport_system_esplicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel,Data_Transport &data_transport, Data_linear_decay &data_linear_decay);


/*!
*Definition of the Transport system solved with an implicit upwind scheme
*/
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel,Data_Transport &data_transport, Data_linear_decay &data_linear_decay);

#endif

