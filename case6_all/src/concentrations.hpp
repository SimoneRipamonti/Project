#ifndef CONCENTRATIONS_HH
#define CONCENTRATIONS_HH

#include <Eigen/Dense>
#include "parameters.hpp"
#include "matrix.hpp"

#include <Eigen/LU>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>

#include <string>


enum Method {ExplicitEuler,PredictorCorrector,Heun}; //Numerical scheme for solving the reaction part

/*!
*\brief Class for the 6 reagents case.
*
*This class has as attributes the concentrations of the reagents and as member functions functions that set and simulate the problem in question.
*/
class Concentration
{
public:
    Concentration(const std::string &filename);/*!<Constructor*/

    unsigned int get_Nx() const;/*!<Getter for Nx*/

    unsigned int get_Nt() const;/*!<Getter for Nt*/

    void set_initial_cond();/*!<Function that sets initial condition for the reagents*/
    /*!
     *Function that sets the matrices for the transport part equation
     *\param M reference to the sparse matrix that solves the transport problem
     *\param rhs reference to the sparse matrix which will be multiplied by the past past solution vector to form the rhs of the transport problem
     *\param vel constant reference to the fluid velocity
    */
    void assemble_transport(Matrix& M, Matrix& rhs, const Vector& vel) const;

    /*!
       *Function that sets the solvers for the transport problems:
       *\param solver reference to the solver used to solve the transport problems for moving reagents (Ca,SiO2,HCO3,H+,CO2)
       *\param solver1 reference to the solver used to solve the transport problem for the immobile solid reagent (CaSiO3)
       *\param M_rhs reference to the sparse mass matrix that will be multiplied by the past past solution vector to form the rhs of the transport problem
       *\param rhs_i reference to the rhs vectors of the i-transport problem
       *\param vel constant reference to the fluid velocity
       *\param Nx spatial steps discretization
    */
    void define_transport_solver(Solver& solver, Solver& solver1, Matrix& M_rhs, Vector& rhs_psi2, Vector& rhs_psi3, const Vector& vel, unsigned int Nx);

    /*!
       *Function that sets the rhs for the total concentrations which have a constant input bc
       *\param rhs_psi2 reference to rhs of the psi_2 total concentration
       *\param rhs_psi3 reference to rhs of the psi_3 total concentration
       *\param vel constant reference to the fluid velocity
       *\param CO2_in input boundary condition for CO2
       *\param CO2_out output boundary condition for CO2
       *\param CO2_bc type of bc
       *\param H_piu_in input boundary condition for H_piu
       *\param H_piu_out output boundary condition for H_piu
       *\param H_piu_bc boundary condition type
    */
    void  assemble_rhs(Vector& rhs_psi2, Vector& rhs_psi3, const Vector& vel, double CO2_in, double CO2_out, const std::string& CO2_bc, double H_piu_in, double H_piu_out, const std::string& H_piu_bc);


    /*!
       *Function that computes the total concentrations from the reagents
       *\param psi_i reference to the i-total concentration at instant equal to t=step*dt
       *\param step is an unsigned int that tells us at what time instant we are t=step*dt
    */
    void compute_psi(unsigned int step, Vector& psi1, Vector& psi2, Vector& psi3, Vector& psi4, Vector& psi5) const;

    /*!
       *Function that computes the reaction term of the equation in the classical standard way with all the reagent parameters (Activation Energy, Area,Temperature, constant reaction term)
       *\param rd reference to the equation reaction term
       *\param step is an unsigned int that tells us at what time instant we are t=step*dt
    */
    void compute_rd(unsigned int step, Vector& rd) const;

    /*!
     * Function that computes the reaction term of the equation using the precipitation constant k_p
     *  \param rd reference to the equation reaction term
     *  \param step is an unsigned int that tells us at what time instant we are t=step*dt
    */
    void compute_rd_kd(unsigned int step, Vector& rd) const;

    /*!
     * Function that computes one temporal step for the transport-reaction equation
     *\param psi_i reference to the i-total concentration
     *\param rd is reference to rhe reaction term of the equation
     *\param step is an unsigned int that tells us at what time instant we are t=step*dt
     *\param rhs reference to the right-hand  mass matrix for the transport equation
     *\param rhs_psi_i reference to to the right-hand vector for the i-
     transport equation
     *\param solver reference to the solver for the transport problem
     *\param solver1 reference to the solver for the CaSiO3, that has not the transport part
    */
    void one_step_transport_reaction(Vector& psi1, Vector& psi2, Vector& psi3, Vector& psi4, Vector& psi5, Vector& rd, const Matrix& rhs, const Vector& rhs_psi2, const Vector& rhs_psi3, unsigned int step, Solver &solver, Solver &solver1);


    /*!
     * Function that solves the reaction part with an Explicit Euler
     *\param psi_i reference to the i-total concentration
     *\param rd is reference to rhe reaction term of the equation
     *\param rhs reference to the right-hand  mass matrix for the transport equation
     *\param rhs_psi_i reference to to the right-hand vector for the i-
     transport equation
     *\param solver reference to the solver for the reagents transport problem 
     *\param solver1 reference to the solver for the CaSiO3, that has not the transport part
    */
    void Euler_Explicit(Vector& psi1, Vector& psi2, Vector& psi3, Vector& psi4, Vector& psi5, const Vector& rd, const Matrix&  rhs, const Vector& rhs_psi2, const Vector& rhs_psi3, Solver &solver, Solver &solver1) const;

    /*!
     * Function that solves the equation for one psi
     *\param psi reference to the total concentration
     *\param M_rhs reference to the right-hand matrix for the transport equation
     *\param rhs_CO2 reference to the vector that identify the rhs for the CO2 transport problem due to the in boundary condition
     *\param rd reference to rhe reaction term of the equation
     *\param solver reference to the solver for the transport problem

    */
    void transport_and_reaction(Vector& psi, const Matrix& M_rhs, const Vector& rhs_CO2, const Vector& rd, Solver &solver) const;

    /*!
     * Functions that solves the non linear system with the Newthon method, in order to compute the real concentrations
     *\param step is an unsigned int that tells us at what time instant we are t=step*dt
     *\param psi_i reference to the i-total concentration
    */

    void compute_concentration(unsigned int step, const Vector& psi1, const Vector& psi2, const Vector& psi3, const Vector& psi4, const Vector& psi5);

    /*!
     * Functions that computes the rhs of the Newton Scheme
     *\param rhs reference to rhe right-hand side to be computed
     *\param old_it reference to the old iterate
     *\param psi_i reference to the i-total concentration
    */
    void compute_rhs(Vector& rhs, const Vector& old_it, double psi1, double psi2, double psi3, double psi4, double psi5) const;

    /*!
     * Functions that computes Jacobian last row
     *\param J reference to the Jacobian Matrix (that is a 6x6 matrix since it is specific for a single spatial node)
     *\param old_it reference to the old iterate
    */
    void compute_Jacob(Matrix_full& J,const Vector& old_it) const;

    /*!
     *Function that saves the solution for a specific reagent in a .csv file in the following way: each row is a spatial position, each column is a time instant
     *\param name string which tells us the reagent
     */
    void output_results_fixed_time(const std::string& name) const;

    /*!
     *Function that saves the solution for a specific reagent in a .csv file in the following way: each row is a time instant, each column is a spatial position
     *\param name string which tells us the reagent
     */
    void output_results_fixed_space(const std::string& name) const;

    /*!
     *Function that saves the solution for each reagent in a specific spatial position during time
    *\param pos is the spatial position where we want to know all the 6 reagents evolution in time
     */
    void output_all_reagents(unsigned int pos) const;

    /*!
     *Function that prints the solution for each reagent at a specific time step
    *\param i temporal step
     */
    void column(unsigned int i);


private:
    Matrix_full Ca;/*!<Matrix that stores the [Ca+] in time and space*/
    Matrix_full H_piu;/*!<Matrix that stores the [H+] in time and space*/
    Matrix_full HCO3_meno;/*!<Matrix that stores the [HCO_3-] in time and space*/
    Matrix_full CO2;/*!<Matrix that stores the [CO_2] in time and space*/
    Matrix_full CaSiO3;/*!<Matrix that stores the [CaSiO_3] in time and space*/
    Matrix_full SiO2;/*!<Matrix that stores the [SiO_2] in time and space*/

    Data_Transport data_transp;/*!<Data for the transport part*/
    Data_6Reagents data_reagents;/*!<Data for the 6 reagents*/
    Data_Reaction data_reaction;/*!<Physical data for the reaction setting*/
    Data_inflow data_CO2;/*!<Input CO2 data (It's the only reagent with a constant inflow*/

    double h;/*!<Spatial step*/
    double dt;/*!<Temporal step*/

    Method method;/*!<Numerical method that solves the reaction part*/
};







#endif
