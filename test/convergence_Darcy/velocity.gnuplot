set datafile separator ','

set key autotitle columnhead # use the first line as title
set ylabel "Pressure" # label for the Y axis
set xlabel 'Space' # label for the X axis



plot "exact_velocity.csv" using 1:2 with lines, "velocity.csv" using 1:2 with lines, 

#plot "exact_velocity.csv" using 1:2 with lines, 
#plot "velocity.csv" using 1:2 with lines, 




