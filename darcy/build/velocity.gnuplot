set term png
set output "velocity.png"

set datafile separator ','

set key autotitle columnhead # use the first line as title
set key font ",10"
set ylabel "Velocity [m/s]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Space [m]' # label for the X axis
set xlabel font ",12"

set yrange [0:0.35]   



plot "exact_velocity.csv" using 1:2 with lines, "velocity10.csv" using 1:2 with lines, "velocity20.csv" using 1:2 with lines, "velocity40.csv" using 1:2 with lines, "velocity80.csv" using 1:2 with lines, "velocity160.csv" using 1:2 with lines,  

#plot "exact_velocity.csv" using 1:2 with lines, 
#plot "velocity.csv" using 1:2 with lines, 



