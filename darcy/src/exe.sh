#!/bin/bash
./main
gnuplot -p velocity.gnuplot 
gnuplot -p pressure.gnuplot 
gnuplot -p error.gnuplot 

