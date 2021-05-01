#ifndef MUPARSER_FUN_HH
#define MUPARSER_FUN_HH

#include <muParser.h>

#include <iostream>
#include <string>

class muparser_fun
{
public:
      muparser_fun(const std::string &);
      double operator() (const double &);
    
private:
    double var;
    mu::Parser parser;
};

#endif

























}
