reset
#set terminal postscript eps enhanced color 14 size 6.5cm,6.5cm

set terminal postscript eps enhanced color 14 size 9cm,6.5cm
set output "latency.eps"
set title "RAM latency (in nanosec)"
# =============================================
set size 1,1
set xlabel ""
set ylabel ""
set yrange [0:70]
set key outside Left reverse
   # set xtics rotate by -45
#set ytics 20
#set key outside right box
#set key outside bottom box
# =============================================
set style data histogram
set style fill pattern border
set style histogram 
# clustered
set grid


plot 'latency.dat' using 4:xticlabels(1) title 'nehalem local' lt -1,\
           '' using 5:xticlabels(1) title 'nehalem remote' lt -1,\
           '' using 6:xticlabels(1) title 'haswell local' lt -1,\
           '' using 7:xticlabels(1) title 'haswell node 1' lt -1,\
           '' using 8:xticlabels(1) title 'haswell node 2' lt -1,\
           '' using 9:xticlabels(1) title 'haswell node 3' lt -1,\
           '' using 2:xticlabels(1) title 'KNL mcdram' lt -1,\
           '' using 3:xticlabels(1) title 'KNL ddr4' lt -1

	       
