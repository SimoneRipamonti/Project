set term png
set output "solutions/CaSiO3.png"

set datafile separator ','

set title 'CaSiO_3'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "[CaSiO_3] [mol/L]" # label for the Y axis
set xlabel 'Space [m]' # label for the X axis
    
plot "CaSiO3_fixed_time.csv" using 1:2 with lines lw 2,'' using 1:4 with lines lw 2,'' using 1:6 with lines lw 2,'' using 1:8 with lines lw 2,'' using 1:10 with lines lw 2,'' using 1:12 with lines lw 2, 






