#!/bin/bash


gnuplot -p -e 'plot "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:5 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:10 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:15 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:20 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:25 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:30 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:35 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:40 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:45 with linespoints, 
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:50 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:55 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:60 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:65 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:70 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:75 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:80 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:85 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:90 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/transport.csv" using 1:95 with linespoints'
