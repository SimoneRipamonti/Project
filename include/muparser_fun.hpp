#ifndef MUPARSER_FUN_HH
#define MUPARSER_FUN_HH

#include <muParser.h>

#include <iostream>
#include <string>


/*!
  *\brief Class for creating and set muparser_fun
  *
  *Class that permit to store functions, passed from the pot file, in muparser function (handle function in Matlab)
*/
class muparser_fun
{
public:
    muparser_fun()=default;/*!<default constructor*/
    muparser_fun(const muparser_fun &m);/*!<copy-constructor*/
    void set_value(const std::string &);/*!<function that takes in input the string that represents the function to be treat and gives as output the respective muparser function*/
    double operator() (const double &);/*!<operator that permit to evaluate the function at a particular value of the independent variable*/
    muparser_fun& operator= (const muparser_fun &p);/*!<copy operator*/
private:
    double var;/*!<Indipendent variable of the muparser function*/
    mu::Parser parser;/*!< Parser from the muparser library*/
};

#endif


























