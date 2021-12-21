set term png
set output "solutions/pressure.png"


set datafile separator ','

set key autotitle columnhead # use the first line as title
set key font ",10"
set ylabel "Pressure [Pa]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Space [m]' # label for the X axis
set xlabel font ",12"



plot "exact_pressure.csv" using 1:2 with lines lw 2, "pressure10.csv" using 1:2 with lines lw 2, "pressure20.csv" using 1:2 with lines lw 2, "pressure40.csv" using 1:2 with lines lw 2, "pressure80.csv" using 1:2 with lines lw 2, "pressure160.csv" using 1:2 with lines lw 2

#plot "exact.csv" using 1:2 with lines, 

#plot "pressure100.csv" using 1:2 with lines, "pressure1600.csv" using 1:2 with lines

#plot "exact.csv" using 1:2 with lines, "pressure100.csv" using 1:2 with lines, "pressure1600.csv" using 1:2 with lines


