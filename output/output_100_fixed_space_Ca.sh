#!/bin/bash


gnuplot -p -e 'plot "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:2 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:11 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:21 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:31 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:41 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:51 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:61 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:71 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:81 with linespoints,
                    "/home/simoripa96/Scrivania/my_project/src/Ca_fixed_space.csv" using 1:91 with linespoints'
