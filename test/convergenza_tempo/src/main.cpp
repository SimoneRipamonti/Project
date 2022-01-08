
#include <cmath>
#include <string>
#include "function.hpp"

#include <iostream>

#include <vector>



int main()
{

    std::vector<unsigned int> Nt{10,20,40,80};
    double T_end{1.};

    double y0{1.0};

    std::vector<double> err_EE(4);
    std::vector<double> err_PC(4);
    std::vector<double> err_IE(4);
    std::vector<double> err_H(4);

    auto source = [](const double &y, const double &t) -> double
    {
        return std::sin(t)*std::pow(y,2);
    };

    auto y_ex = [y0](const double &t) -> double
    {
        return -y0/(y0-y0*std::cos(t)-1);
    };

    for (unsigned int i=0; i<Nt.size(); ++i)

    {
        double dt{static_cast<double>(T_end/Nt[i])};

        std::vector<double> y_EE(Nt[i]+1, 0.0);
        std::vector<double> y_PC(Nt[i]+1, 0.0);
        std::vector<double> y_IE(Nt[i]+1, 0.0);
        std::vector<double> y_H(Nt[i]+1, 0.0);

        y_EE[0]=y0;
        y_PC[0]=y0;
        y_IE[0]=y0;
        y_H[0]=y0;

        for(unsigned int tt=0; tt<Nt[i]; ++tt)
        {
            y_EE[tt+1]=EE(tt,dt,y_EE[tt],source);
            y_PC[tt+1]=predictor_corrector(tt,dt,y_PC[tt],source);
            y_IE[tt+1]=IE(tt,dt,y_IE[tt],source);
            y_H[tt+1]=Heun(tt,dt,y_H[tt],source);
        }


        err_EE[i]=std::abs(y_EE.back()-y_ex(T_end));
        err_PC[i]=std::abs(y_PC.back()-y_ex(T_end));
        err_IE[i]=std::abs(y_IE.back()-y_ex(T_end));
        err_H[i]=std::abs(y_H.back()-y_ex(T_end));

    }


    double p_EE {std::log(err_EE[2]/err_EE[3])/std::log(2)};
    double p_PC {std::log(err_PC[2]/err_PC[3])/std::log(2)};
    double p_IE {std::log(err_IE[2]/err_IE[3])/std::log(2)};
    double p_H  {std::log(err_H[2]/err_H[3])/std::log(2)};

    std::cout<<"Convergence order for EE: "<<p_EE<<std::endl;
    std::cout<<"Convergence order for PC: "<<p_PC<<std::endl;
    std::cout<<"Convergence order for IE: "<<p_IE<<std::endl;
    std::cout<<"Convergence order for H: "<<p_H<<std::endl;
    


}



