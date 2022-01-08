This folder contains the header files of functions and classes which are used in the different case folders.

We define 6 modules (groups) in order to document better (with doxygen) the content of these files.

# Matrices
It groups all the classes defined in the matrix.hpp file. The matrices that we are talking about are the ones obtained by the discretization of the equations.

First of all abstract class (the Abstract Matrix class) it's declared which is the prototype for all the derived class whcich inherit from it. 

The derived classes are:

Matrix A is the mass velocity matrix for the Darcy system;

Matrix B is the saddle point matrix of the Darcy system;

Matrix C is the mass matrix obtained in the transport equation with the finite volume method; 

Matrix F_piu_ and F_meno are the ones that define the upwind scheme;

Matrix R is the reactive matrix used in the case5_2reagents for defining the reactive part of the equation.

# Parameters
It groups all the classes defined in the parameters.hpp file. 

These classes are used to initialize all the data that are used in order to define the different problems.

Their constructors take as input a data.pot file and through GetPot they store the data in different variables. 

# Darcy Functions
It groups all the functions in the darcy.hpp file. 

These functions are the ones that are used to define and solve the Darcy Problem.

# Darcy Output
It groups all the functions in the darcy_output.hpp file. 

These functions print on csv files the solution given by the code and the exact ones given by a suitable muparser function.

There is also a function that is used to print the error of the convergence analysis in order to see the convergence order of the scheme adopted. 

# Output Function
It groups all the function in the output.hpp file. 

These functions print on csv files the solution of temporal and spatial problems such as transport and reaction ones. 

# MuParser
It contains the muparser_fun_ class which is defined in the muparser_fun.hpp file.

This class permit us to read and use functions which are given as handle function in the data.pot file. 

# Functions
It contains 3 different functions:

1. the first one is specific for setting the solver to a linear system

2. the second and the third ones are used in the transport_decay.cpp file in order to set the linear system for the transport and decay problem

# Transport and decay functions
It contains the functions which define and solve the transport and decay problem.



















