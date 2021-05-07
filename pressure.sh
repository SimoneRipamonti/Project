#!/bin/bash

gnuplot -p -e 'plot "pressure.csv" using 1:2 with linespoints'
