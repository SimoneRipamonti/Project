#!/bin/bash


gnuplot -p -e 'plot "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:2 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:3 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:4 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:5 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:6 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:7 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:8 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:9 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:10 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/transport.csv" using 1:11 with linespoints'
