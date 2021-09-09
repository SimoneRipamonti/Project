set datafile separator ','

set title 'Ca'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "[Ca^{2+}] [mol/L]" # label for the Y axis
set xlabel 'Time [s]' # label for the X axis

#set yrange [0:1.5]

set logscale x
set logscale y                    
  
plot "Ca_fixed_space.csv" using 1:2 with lines,'' using 1:4 with lines,'' using 1:6 with lines,'' using 1:8 with lines,'' using 1:10 with lines,'' using 1:12 with lines, 






