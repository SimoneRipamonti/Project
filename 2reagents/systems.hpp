#ifndef SYSTEMS_HH
#define SYSTEMS_HH

#include <Eigen/Dense>
#include "matrix.hpp"
#include "parameters.hpp"

//Definition of the Darcy Syste
//void set_Darcy_system(Data_Darcy &data, Eigen::MatrixXd &Matrix,Eigen::VectorXd &rhs);

/*!
 *Definition of the Tranport-Reaction equation for 2 reagents solved with the following schemse: explicit for both transport and reaction
 */
void Transport_system_esplicit_2_reagents(Eigen::MatrixXd &Ca,Eigen::MatrixXd &CaSiO3, Eigen::VectorXd &vel,Data_Transport &data_transport,Data_Reaction &data_reaction, Data_2Reagents &data_2reagents);

/*!
 *Definition of the Tranport-Reaction equation for 2 reagents solved with the following schemse: implicit for transport and explicit for reaction
 */
void Transport_system_implicit_2_reagents(Eigen::MatrixXd &Ca,Eigen::MatrixXd &CaSiO3,Eigen::VectorXd &velocity,Data_Transport &data_transport, Data_Reaction &data_reaction, Data_2Reagents &data_2reagents);

#endif

