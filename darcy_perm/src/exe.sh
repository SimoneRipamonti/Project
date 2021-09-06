#!/bin/bash

./main
gnuplot -p velocity.gnuplot 
gnuplot -p pressure.gnuplot

