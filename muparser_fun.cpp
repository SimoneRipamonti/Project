#include "muparser_fun.hpp"

muparser_fun::muparser_fun(const std::string &s):p(m.p)
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
       };
double muparser_fun::operator()(const double &x)
{
      double y;
      var = x;
      try
        {
          y=p.Eval();
        }    
      catch(mu::Parser::exception_type &e)
        {
          std::cerr<<e.GetMsg()<<std::endl;
        }
      return y;
};
     



