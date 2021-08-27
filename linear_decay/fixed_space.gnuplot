set datafile separator ','

set title 'FIXED SPACE'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "Concentration C" # label for the Y axis
set xlabel 'Time' # label for the X axis


#plot "Ca_fixed_space.csv" using 1:2 with lines,'' using 1:3 with lines,'' using 1:4 with lines,'' using 1:5 with lines, '' using 1:6 with lines, '' using 1:7 with lines, '' using 1:8 with lines, '' using #1:9 with lines, '' using 1:10 with lines, 

plot "Ca_fixed_space.csv" using 1:2 with lines,'' using 1:4 with lines,'' using 1:6 with lines,'' using 1:8 with lines,'' using 1:10 with lines, 


