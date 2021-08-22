set datafile separator ','

set key autotitle columnhead # use the first line as title
set ylabel "Pressure" # label for the Y axis
set xlabel 'Space' # label for the X axis



plot "exact_pressure.csv" using 1:2 with lines, "pressure.csv" using 1:2 with lines

#plot "exact.csv" using 1:2 with lines, 

#plot "pressure100.csv" using 1:2 with lines, "pressure1600.csv" using 1:2 with lines

#plot "exact.csv" using 1:2 with lines, "pressure100.csv" using 1:2 with lines, "pressure1600.csv" using 1:2 with lines


