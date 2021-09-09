#!/bin/bash

./main
mkdir -p solutions
gnuplot -p fixed_time.gnuplot 
gnuplot -p fixed_space.gnuplot

