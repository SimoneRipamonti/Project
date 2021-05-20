#!/bin/bash

gnuplot -p -e 'plot "/home/simoripa96/Scrivania/my_project/pressure.csv" using 1:2 with linespoints'
