set datafile separator ','

set title 'SiO_2'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "[SiO_2] [mol/L]" # label for the Y axis
set xlabel 'Time [s]' # label for the X axis

#set yrange [0:1.5]                     # intervallo x in questo caso 0-5

set logscale x
set logscale y    


plot "SiO2_fixed_space.csv" using 1:2 with lines,'' using 1:4 with lines,'' using 1:6 with lines,'' using 1:8 with lines,'' using 1:10 with lines,'' using 1:12 with lines, 






