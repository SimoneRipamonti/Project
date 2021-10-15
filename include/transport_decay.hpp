#ifndef TRANSPORT_DECAY_HH
#define TRANSPORT_DECAY_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"


/** \addtogroup Transport_decay_functions
 *\brief Functions that implement and solve the transport and decay system
 *  @{
*/
/*!
 *Function that implements and solves the explicit transport and decay system (decay part is treated explicitly)
 *\param Ca is the matrix solution where we store our solution at each istant: each row represent a spatial position, each column represent a time istant.
 *\param vel is the velocity vector which can be computed by the Darcy problem
 *\param data_transport stores the data related to the transport part
 *\param data_linear_decay stores the data related to the linear decay part
*/
void Transport_system_explicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &vel, Data_Transport &data_transport, Data_linear_decay &data_linear_decay);

/*!
 *Function that implements and solves the implicit transport and decay system (decay part is treated explicitly)
 *\param Ca is the matrix solution where we store our solution at each istant: each row represent a spatial position, each column represent a time istant.
 *\param vel is the velocity vector which can be computed by the Darcy problem
 *\param data_transport stores the data related to the transport part
 *\param data_linear_decay stores the data related to the linear decay part
*/
void Transport_system_implicit(Eigen::MatrixXd &Ca, Eigen::VectorXd &velocity,Data_Transport &data_transport, Data_linear_decay &data_linear_decay);

/** @}*/
#endif

