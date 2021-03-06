#include <iostream>
#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"
#include "output_darcy.hpp"
#include "darcy.hpp"
#include "functions.hpp"
#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <cmath>



int main(int argc, char **argv)
{
    Data_Darcy data("data.pot");

    Eigen::VectorXi a(5); //vector that contains the steps number
    Vector errp(5);//vector pressure error
    Vector errvel(5);//vector velocity error

    a<<10,20,40,80,160;


    muparser_fun sol_ex; //muparser_fun that contains the exact solution for pressure
    sol_ex.set_value("1e7*(1/(_pi*_pi*16)*sin(_pi*4*x)-0.1*x-1/(4*_pi)*x+0.1)");

    muparser_fun sol_vel; //muparser_fun that contains the exact solution for velocity
    sol_vel.set_value("-1/(4*_pi)*cos(4*_pi*x)+0.1+1/(4*_pi)");


    Vector sol;//vector that stores the numerical solution

    Vector exact; //vector that stores the exact solution for pressure
    Vector exact_vel; //vector that stores the exact solution for velocity
   
    //Loop where at each iterate we change the number of steps
    for (unsigned int i=1; i<6; ++i)
    {
        data.Nx=a(i-1);

        double h =static_cast<double>(data.L)/data.Nx; //space step

        sol.resize(data.Nx+data.Nx+1); //solution vectors are resized
        exact.resize(data.Nx);
        exact_vel.resize(data.Nx+1);

        Matrix M(data.Nx+data.Nx+1,data.Nx+data.Nx+1);//Initialization of the big matrix for the Darcy system
        Vector rhs(data.Nx+data.Nx+1);//Initialization of the rhs of Darcy
        set_Darcy_system(data,M,rhs,h);//Definition of the Darcy system Mx=rhs

        Solver solver; //solver for the linear system
        set_solver(M,solver);

        sol= solver.solve(rhs);//The Darcy system is solved and the solution is stored in the sol vector

       //The exact solution vectors are filled
        for(unsigned int j=0; j<data.Nx; ++j)
            exact(j)=sol_ex(h/2+j*h);
        for(unsigned int j=0; j<data.Nx+1; ++j)
            exact_vel(j)=sol_vel(j*h);

        Darcy_output_results(sol,data.Nx,data.L);

        //Errors computation
        errp(i-1)=(exact-sol.tail(data.Nx)).norm();
        errvel(i-1)=(exact_vel-sol.head(data.Nx+1)).norm();

    }

    std::cout<<"Order of convergence for pressure:"<<std::log(errp(3)/errp(4))/std::log(2.0)<<std::endl;
    
    std::cout<<"Order of convergence for velocity:"<<std::log(errvel(3)/errvel(4))/std::log(2.0)<<std::endl;


//Plot of the output results
    Darcy_output_results(sol,data.Nx,data.L);
    pressure_exact_result(exact,data.Nx,data.L);
    velocity_exact_result(exact_vel,data.Nx,data.L);
    output_error(errp,a,"pressure");
    output_error(errvel,a,"velocity");




    return 0;
}
