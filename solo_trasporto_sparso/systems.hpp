#ifndef SYSTEMS_HH
#define SYSTEMS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"



//Definition of the Tranport system solved with an esplicit scheme
void Transport_system_esplicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport);

//Definition of the Tranport system solved with an implicit scheme
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &velocity,Data_Transport &data_transport);

#endif

