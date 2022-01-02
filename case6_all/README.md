#Case6: Final overall case with 6 reagents
This test-case implements and solve the overall case introducing a new class which is the Concentration class.

The reactions involved in the process simulated are:

1. The acid carbonic dissolution:

        H2CO3 <------> H+ + HCO_3-

which is treated as an equilibrium reaction since it's much faster then the CaSiO3 dissolution. 
 
2. Wollastonite dissolution:

         CaSiO3+2H+  --------->  Ca2+ + SiO2(aq) + H_2O

which is treated as a cinetic reaction with its own reaction rate rd.
