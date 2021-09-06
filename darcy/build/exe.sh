#!/bin/bash
./my_exec
gnuplot -p velocity.gnuplot 
gnuplot -p pressure.gnuplot 
gnuplot -p error.gnuplot 

