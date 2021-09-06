set datafile separator ','

set title 'CaSiO_3'
set title font " ,12"
set key autotitle columnhead # use the first line as title
set ylabel "CaSiO_3[mol/L]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Time[s]' # label for the X axis
set xlabel font ",12"



#plot "CaSiO3_fixed_space.csv" using 1:2 with lines,'' using 1:4 with lines,'' using 1:6 with lines,'' using 1:8 with lines,'' using 1:10 with lines, 
plot "CaSiO3_fixed_space.csv"  using 1:12 with lines,






