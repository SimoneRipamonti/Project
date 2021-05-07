#!/bin/bash

gnuplot -p -e 'plot "velocity.csv" using 1:2 with linespoints'
