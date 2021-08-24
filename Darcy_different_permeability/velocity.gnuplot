set datafile separator ','

set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "Velocity" # label for the Y axis
set xlabel 'Space' # label for the X axis



#plot "exact_velocity.csv" using 1:2 with lines, "velocity10.csv" using 1:2 with lines, "velocity20.csv" using 1:2 with lines, "velocity40.csv" using 1:2 with lines, "velocity80.csv" using 1:2 with lines, #"velocity160.csv" using 1:2 with lines,  

#plot "exact_velocity.csv" using 1:2 with lines, 
plot "velocity100.csv" using 1:2 with lines, 




