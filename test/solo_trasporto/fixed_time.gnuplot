set datafile separator ','

set key autotitle columnhead # use the first line as title
set ylabel "First Y Units" # label for the Y axis
set xlabel 'Space' # label for the X axis

#set xrange [0:1]                     # intervallo x in questo caso 0-5
  




#plot "Ca_fixed_time.csv" using 1:2 with lines,'' using 1:11 with lines,'' using 1:21 with lines,'' #using 1:31 with lines, '' using 1:41 with lines, '' using 1:51 with lines, '' using 1:61 with #lines, #'' using 1:71 with lines, '' using 1:81 with lines, '' using 1:91 with lines

plot "Ca_fixed_time.csv" using 1:2 with lines,'' using 1:3 with lines,'' using 1:4 with lines,'' using 1:5 with lines, '' using 1:6 with lines, '' using 1:7 with lines, '' using 1:8 with lines, '' using 1:9 with lines, '' using 1:10 with lines, 

#plot "Ca_fixed_time.csv" using 1:2 with lines,''




