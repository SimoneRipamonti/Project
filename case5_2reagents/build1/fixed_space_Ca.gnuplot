set term png
set output "solutions/fixed_space_Ca.png"

set datafile separator ','

set title 'Ca^{2+}'
set title font " ,12"
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "Ca^{2+}[mol/L]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Time[s]' # label for the X axis
set xlabel font ",12"


#plot "Ca_fixed_space.csv" using 1:2 with lines,'' using 1:3 with lines,'' using 1:4 with lines,'' using 1:5 with lines, '' using 1:6 with lines, '' using 1:7 with lines, '' using 1:8 with lines, '' using 1:9 with lines, '' using 1:10 with lines, 

#plot "Ca_fixed_space.csv" using 1:2 with lines,'' using 1:4 with lines,'' using 1:6 with lines,'' using 1:8 with lines,'' using 1:10 with lines,'' using 1:12 with lines, 

plot "Ca_fixed_space.csv" using 1:12 with lines,



