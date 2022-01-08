set term png
set output "solutions/fixed_space.png"

set datafile separator ','

set title 'FIXED SPACE'
set title font ",12"
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "Concentration [mol/L]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Time [s]' # label for the X axis
set xlabel font ",12"

set yrange [0:1.5]

plot "Ca_fixed_space.csv" using 1:2 with lines lw 2,'' using 1:4 with lines lw 2,'' using 1:6 with lines lw 2,'' using 1:8 with lines lw 2,'' using 1:10 with lines lw 2,'' using 1:12 with lines lw 2,



