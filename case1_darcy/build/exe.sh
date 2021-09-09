#!/bin/bash

./main

mkdir -p solutions
gnuplot -p velocity.gnuplot 
gnuplot -p pressure.gnuplot 
gnuplot -p error.gnuplot 

