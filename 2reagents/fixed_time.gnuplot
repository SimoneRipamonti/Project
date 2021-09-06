set datafile separator ','

set title 'Ca^{2+}'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "[Ca^{2+}] [mol/L]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Space [m]' # label for the X axis
set xlabel font ",12"

                     # intervallo x in questo caso 0-5
  
plot "Ca_fixed_time.csv" using 1:2 with lines,'' using 1:4 with lines,'' using 1:6 with lines,'' using 1:8 with lines,'' using 1:10 with lines,'' using 1:12 with lines, 






