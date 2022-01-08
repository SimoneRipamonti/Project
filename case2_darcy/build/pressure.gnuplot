set term png
set output "solutions/pressure.png"

set datafile separator ','

set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "Pressure [Pa]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Space [m]' # label for the X axis
set xlabel font ",12"

plot "pressure100.csv" using 1:2 with lines lw 2,



    
