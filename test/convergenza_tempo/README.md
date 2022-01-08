This folder presents a convergence order analysis for 4 different numerical schemes used to solve an ODE. The program prints as output te convergence of order for the different test cases.

The example problem tested is the following one:
	
	dy/dt=f(t,y) with f(t,y)=sin(t)y^2 and y(0)=y0

This differential problem has as exact solution:

	y_ex=-y0/(y0-y0cos(t)-1)

The 4 order schemes tested are the following ones:
1. Explicit Euler (EE)
2. Implicit Euler (IE)
3. Predictor-Corrector (PC)
4. Heun (HE)

