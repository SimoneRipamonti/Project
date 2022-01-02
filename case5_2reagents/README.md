#Case5: 2 reagents
This test case implements and solves two coupled PDEs, which coupling is due to the reaction term. 

In particular the reaction involved is CaSiO_3 dissolution:

CaSiO_3 + 2H^+ ---->Ca^2+ + SiO_2(aq) + H_2O

Here they are implemented:

1. The transport and reaction equation for the $Ca^{2+}$ ion (a PDE)
2. The reaction equation for $CaSiO_3$ (an ODE, there is not the transport part here since the mineral is immobile)
