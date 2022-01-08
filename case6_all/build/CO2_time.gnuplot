set term png
set output "solutions/CO2.png"

set datafile separator ','

set title 'CO_2'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "CO_2 [mol/L]" # label for the Y axis
set xlabel 'Space [m]' # label for the X axis
   
plot "CO2_fixed_time.csv" using 1:2 with lines lw 2,'' using 1:4 with lines lw 2,'' using 1:6 with lines lw 2,'' using 1:8 with lines lw 2,'' using 1:10 with lines lw 2,'' using 1:12 with lines lw 2






