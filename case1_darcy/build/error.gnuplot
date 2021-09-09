set term png
set output "solutions/error.png"

set datafile separator ','

set key autotitle columnhead # use the first line as title
set ylabel "Error" # label for the Y axis
set xlabel 'Nx' # label for the X axis

set logscale x
set logscale y

plot "error.csv" using 1:2 with lines, "error.csv" using 1:3 with lines,

