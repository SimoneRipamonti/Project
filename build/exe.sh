#!/bin/bash

./main
gnuplot -p fixed_time.gnuplot 
gnuplot -p fixed_space.gnuplot

