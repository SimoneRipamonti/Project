set datafile separator ','

set key autotitle columnhead # use the first line as title
set ylabel "First Y Units" # label for the Y axis
set xlabel 'Space' # label for the X axis


plot "Ca_fixed_time.csv" using 1:2 with lines,'' using 1:11 with lines,'' using 1:21 with lines,'' using 1:31 with lines, '' using 1:41 with lines, '' using 1:51 with lines, '' using 1:61 with lines, '' using 1:71 with lines, '' using 1:81 with lines, '' using 1:91 with lines





