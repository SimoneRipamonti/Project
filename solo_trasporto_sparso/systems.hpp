#ifndef SYSTEMS_HH
#define SYSTEMS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"


//void Transport_system_esplicit(Eigen::MatrixXd &solution,Eigen::VectorXd &velocity,Data_Transport &data);


void Transport_system_esplicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport, Data_2Reagents &data_2reagents);

//Definition of the Tranport system solved with an implicit scheme
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &velocity,Data_Transport &data_transport, Data_2Reagents &data_2reagents);

#endif

