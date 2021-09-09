#ifndef SYSTEMS_HH
#define SYSTEMS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"


/*!
 *Definition of the Tranport system solved with an explicit scheme
 */
void Transport_system_esplicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport, Data_linear_decay &initial_cond);

/*!
 *Definition of the Tranport system solved with an implicit scheme
*/
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &velocity,Data_Transport &data_transport, Data_linear_decay &initial_cond);

#endif

