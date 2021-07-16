#ifndef MUPARSER_FUN_HH
#define MUPARSER_FUN_HH

#include <muParser.h>

#include <iostream>
#include <string>

//Class that permit to store function, passed from the pot file, in muparser function (handle function in Matlab)
class muparser_fun
{
public:
    muparser_fun()=default;//default constructor
    muparser_fun(const muparser_fun &m);//copy-constructor
    void set_value(const std::string &);//set value of the parser member
    double operator() (const double &);//operator that permit to evaluate the function at a particular value of the independent variable
    muparser_fun& operator= (const muparser_fun &p);//equal operator
private:
    double var;
    mu::Parser parser;
};

#endif


























