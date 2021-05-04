#ifndef MUPARSER_FUN_HH
#define MUPARSER_FUN_HH

#include <muParser.h>

#include <iostream>
#include <string>

class muparser_fun
{
public:
      muparser_fun()=default;
      muparser_fun(const muparser_fun &m);
      void set_value(const std::string &);
      double operator() (const double &);
      muparser_fun& operator= (const muparser_fun &p);
private:
    double var;
    mu::Parser parser;
};

#endif


























