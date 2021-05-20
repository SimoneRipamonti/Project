#!/bin/bash

gnuplot -p -e 'plot "transport.csv" using 1:2 with linespoints,
                    "transport.csv" using 1:6 with linespoints,
                    "transport.csv" using 1:11 with linespoints'
