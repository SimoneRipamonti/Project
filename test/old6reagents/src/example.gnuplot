set datafile separator ','

set key autotitle columnhead # use the first line as title
set ylabel "First Y Units" # label for the Y axis
set xlabel 'Time' # label for the X axis

set logscale x
set logscale y

plot "all.csv" using 1:2 with lines,'' using 1:3 with lines,'' using 1:4 with lines,'' using 1:5 with lines


