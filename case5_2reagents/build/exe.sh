#!/bin/bash

./main
mkdir -p solutions
gnuplot -p fixed_time_Ca.gnuplot 
gnuplot -p fixed_time_CaSiO3.gnuplot
gnuplot -p fixed_space_Ca.gnuplot 
gnuplot -p fixed_space_CaSiO3.gnuplot 
