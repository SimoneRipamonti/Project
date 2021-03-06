#include "muparser_fun.hpp"

muparser_fun::muparser_fun(const muparser_fun &m)//Copy-constructor definition
    : parser(m.parser)
{
    parser.DefineVar("x", &var);
}


void muparser_fun::set_value(const std::string &s)//Definition of the set_value() function
{
    try
    {
        parser.DefineVar("x",&var);
        parser.SetExpr(s);
    }
    catch(mu::Parser::exception_type &e)
    {
        std::cerr<<e.GetMsg()<<std::endl;
    }
}

double muparser_fun::operator()(const double &x)//Definition of the () operator
{
    double y;
    var = x;
    try
    {
        y=parser.Eval();
    }
    catch(mu::Parser::exception_type &e)
    {
        y=0;
        std::cerr<<e.GetMsg()<<std::endl;
    }
    return y;
}

muparser_fun& muparser_fun::operator= (const muparser_fun &p)//Defintion of the = operator
{
    if(&p!=this)
    {

        parser=p.parser;
        parser.DefineVar("x",&var);
    }
    return *this;
}



