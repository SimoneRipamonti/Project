#!/bin/bash

gnuplot -p -e 'plot "output.csv" using 1:2 with linespoints,
                    "output.csv" using 1:3 with linespoints,
                    "output.csv" using 1:4 with linespoints,
		    "output.csv" using 1:5 with linespoints,
		    "output.csv" using 1:6 with linespoints,
		    "output.csv" using 1:7 with linespoints,
	            "output.csv" using 1:8 with linespoints,
		    "output.csv" using 1:9 with linespoints,
	            "output.csv" using 1:10 with linespoints,
                    "output.csv" using 1:11 with linespoints'
