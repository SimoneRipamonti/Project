set term png
set output "solutions/H_piu.png"

set datafile separator ','

set title 'H^+'
set key autotitle columnhead # use the first line as title
set key font ",12"
set ylabel "[H^+] [mol/L]" # label for the Y axis
set xlabel 'Space [m]' # label for the X axis


#set yrange [0:1.5]                     # intervallo x in questo caso 0-5

#set logscale x
#set logscale y    
  
plot "H_piu_fixed_time.csv" using 1:2 with lines lw 2,'' using 1:4 with lines lw 2,'' using 1:6 with lines lw 2,'' using 1:8 with lines lw 2,'' using 1:10 with lines lw 2,'' using 1:12 with lines lw 2, 






