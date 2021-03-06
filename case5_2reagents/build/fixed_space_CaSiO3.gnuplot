set term png
set output "solutions/fixed_space_CaSiO3.png"

set datafile separator ','

set title 'CaSiO_3'
set title font " ,12"
set key autotitle columnhead # use the first line as title
set ylabel "CaSiO_3[mol/L]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Time[s]' # label for the X axis
set xlabel font ",12"

plot "CaSiO3_fixed_space.csv"  using 1:12 with lines lw 2,






