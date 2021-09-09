#!/bin/bash

./main
mkdir -p solutions
gnuplot -p Ca_time.gnuplot 
gnuplot -p CaSiO3_time.gnuplot 
gnuplot -p CO2_time.gnuplot 
gnuplot -p H_piu_time.gnuplot 
gnuplot -p HCO3_meno_time.gnuplot 
gnuplot -p SiO2_time.gnuplot
