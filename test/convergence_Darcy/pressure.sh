#!/bin/bash

gnuplot -p -e 'plot "/home/simoripa96/Scrivania/my_project/convergence/pressure.csv" using 1:2 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/convergence/exact.csv" using 1:2 with linespoints'