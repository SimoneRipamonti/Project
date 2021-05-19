#!/bin/bash

gnuplot -p -e 'plot "output.csv" using 1:2 with linespoints'
