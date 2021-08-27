set datafile separator ','

set key autotitle columnhead # use the first line as title
set ylabel "First Y Units" # label for the Y axis
set xlabel 'Time' # label for the X axis

set logscale x
set logscale y

plot "CaSiO3_fixed_space.csv" using 1:2 with lines


