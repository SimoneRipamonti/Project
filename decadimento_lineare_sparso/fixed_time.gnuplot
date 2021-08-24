set datafile separator ','

set title 'FIXED TIME'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "Concentration C" # label for the Y axis
set xlabel 'Space' # label for the X axis

set yrange [0:1.5]                     # intervallo x in questo caso 0-5
  
plot "Ca_fixed_time.csv" using 1:2 with lines,'' using 1:4 with lines,'' using 1:6 with lines,'' using 1:8 with lines,'' using 1:10 with lines, 






