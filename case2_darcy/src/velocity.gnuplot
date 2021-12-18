set term png
set output "solutions/velocity.png"

set datafile separator ','

set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "Velocity [m/s]" # label for the Y axis
set ylabel font ",12"
set xlabel 'Space [m]' # label for the X axis
set xlabel font ",12"



##plot "exact_velocity.csv" using 1:2 with lines lw 2, "velocity10.csv" using 1:2 with lines lw 2, "velocity20.csv" using 1:2 with lines lw 2, "velocity40.csv" using 1:2 with lines lw 2, "velocity80.csv" using 1:2 with lines lw 2, #"velocity160.csv" using 1:2 with lines lw 2,  

#plot "exact_velocity.csv" using 1:2 with lines, 
plot "velocity100.csv" using 1:2 with lines lw 2, 




