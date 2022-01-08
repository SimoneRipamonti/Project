#include "function.hpp"



double EE (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f )
{
    double y_p{reaction_rate(tt,dt,y_old,f)};
    return y_old+dt*y_p;
}

double predictor_corrector (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f)
{
    double y_p{reaction_rate(tt,dt,y_old,f)};
    double y_star{y_old+dt*y_p};
    double y_p2{reaction_rate(tt+1,dt,y_star,f)};
    return y_old+dt*y_p2;

}

double Heun (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f)
{
    double y_p{reaction_rate(tt,dt,y_old,f)};
    double y_star{y_old+dt*y_p};
    double y_p2{reaction_rate(tt+1,dt,y_star,f)};
    return y_old+0.5*dt*(y_p+y_p2);
}

double IE (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f)
{
    double tol{1.0e-10};
    double err{tol+1};
    unsigned int maxit{1000};
    unsigned int kk{0};
    double y_k{y_old};

    double y_p{0.};
    double J{0.};
    double res{0.};
    double dy{0.};

    for(; kk<maxit and err>tol; ++kk)
    {
        y_p=reaction_rate(tt,dt,y_k,f);
        res=y_k-y_old-dt*y_p;
        J=compute_J(dt,tt,y_k);
        dy=-res/J;
        err=std::abs(dy);
        y_k+=dy;
    }

    if (kk==maxit and err>tol)
        std::cout<<"Mancata convergenza"<<std::endl;

    return y_k;
}


double compute_J(double dt, unsigned int tt, double y)
{
    return 1-dt*2*y*std::sin(tt*dt);

}

double reaction_rate (unsigned int tt, double dt, double y_old, const std::function<double(const double &, const double &)> &f)
{
    return f(y_old,tt*dt);
}


void output_results(const std::vector<unsigned int>& time, const std::vector<double>& y_EE, const std::vector<double>& y_IE, const std::vector<double>& y_PC, const std::vector<double>& y_H,const std::vector<double>& y_ex)
{
    // Output results to CSV file.
    std::ofstream file("output.csv", std::ofstream::out);
    file << "time, EE, IE, PC, H, EX" << std::endl;

    for (unsigned int step = 0; step <time.size(); ++step)
    {
        file << time[step] << ", " << y_EE[step] << ", "
             << y_IE[step] << ", " << y_PC[step]<< ", "<<y_H[step]<<", "<<y_ex[step]<<std::endl;
    }
    file.close();

    // Plot results.
    Gnuplot gp;
    gp << "set xlabel 'Time [days]'; set ylabel 'No. of people'; set key center "
       "right; plot "
       << gp.file1d(std::tie(time, y_EE))
       << "with line linewidth 2 title 'y_EE',"
       << gp.file1d(std::tie(time, y_IE))
       << "with line linewidth 2 title 'y_IE',"
       << gp.file1d(std::tie(time, y_PC))
       << "with line linewidth 2 title 'y_PC',"
       << gp.file1d(std::tie(time, y_H))
       << "with line linewidth 2 title 'y_H',"
       << gp.file1d(std::tie(time, y_ex))
       << "with line linewidth 2 title 'y_Ex',"<< std::endl;
}
