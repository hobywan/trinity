#reset
##set terminal postscript eps enhanced color 14 size 6.5cm,6.5cm
#
#kernel=3
#set terminal postscript eps enhanced color 14 size 8cm,6.5cm
#set output "repart_".kernel.".eps"
#set title "Time ratio per step [Kernel ".kernel."]"
## =============================================
#set size 1,1
#set xlabel ""
#set ylabel ""
#set yrange [0:100]
#   # set xtics rotate by -45
##set ytics 20
##set key outside right box
##set key outside bottom box
## =============================================
#set style data histogram
#unset key
##set style histogram rowstack gap 1
#set style fill solid border -1
#set style fill solid 1
##set boxwidth 0.25 relative
#set boxwidth 0.25
#set style line 2 lc rgb 'black' lt 1 lw 1
#set grid y
#set style histogram rowstacked
#set xtics rotate by 45 right
##set style fill pattern border -1
##iset ytics 10 nomirror
##set logscale y
#file = '../benchs/repart_'.kernel.'.dat' 
#plot newhistogram '',\
#file using 2:xticlabels(1) title 'step 1' fillstyle pattern 1  lt -1,\
#'' using 3:xticlabels(1) title 'step 2' fillstyle pattern 3  lt -1,\
#'' using 4:xticlabels(1) title 'step 3' fillstyle pattern 7  lt -1,\
#'' using 5:xticlabels(1) title 'step 4' fillstyle pattern 10 lt -1  
#
# =============================================
reset
set terminal postscript eps enhanced color 14 size 8cm,6cm
set title "Intel Haswell"
# =============================================
#set size 1,1
set xlabel "threads"
set ylabel "Ratio (%)"
set yrange [0:100]
   # set xtics rotate by -45
#set ytics 20
#set key outside right box
#set key outside bottom box
# =============================================
set style data histogram
#unset key
#set key outside left
set key outside top center
set size ratio 0.5
#set style histogram rowstack gap 1
set style fill solid border -1
set style fill solid 1
#set boxwidth 0.25 relative
set boxwidth 0.25
set style line 2 lc rgb 'black' lt 1 lw 1
set grid y
set style histogram rowstacked

#do for [k=1:4] {
#  file = '_'.k.'/perf_hsw_shock.dat' 
#  set output "hsw_shock_repar_".k.".eps"
#  plot newhistogram '',\
#  file using ($10*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'step 1' fillstyle pattern 1  lt -1,\
#    '' using ($11*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'step 2' fillstyle pattern 3  lt -1,\
#    '' using ($12*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'step 3' fillstyle pattern 7  lt -1,\
#    '' using ($13*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'step 4' fillstyle pattern 10 lt -1,\
#    '' using ($14*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'step 5' fillstyle pattern 11 lt -1  	       
#}

file = '_1/perf_hsw_shock.dat' 
set output "plots/hsw_shock_repar_1.eps"
plot newhistogram '',\
file using ($10*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'filter' fillstyle pattern 1  lt -1,\
  '' using ($11*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'steiner'fillstyle pattern 2  lt -1,\
  '' using ($12*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'kernel' fillstyle pattern 3  lt -1,\
  '' using ($13*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'repair' fillstyle pattern 2 lt rgb "#8B0000"

file = '_2/perf_hsw_shock.dat' 
set output "plots/hsw_shock_repar_2.eps"
plot newhistogram '',\
file using ($11*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'primal' fillstyle pattern 1 lt -1,\
  '' using ($10*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'filter' fillstyle pattern 2 lt -1,\
  '' using ($12*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'indep'  fillstyle pattern 3 lt rgb "#8B0000",\
  '' using ($13*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'kernel' fillstyle pattern 3 lt -1,\
  '' using ($14*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'repair' fillstyle pattern 2 lt rgb "#8B0000"

file = '_3/perf_hsw_shock.dat' 
set output "plots/hsw_shock_repar_3.eps"
plot newhistogram '',\
file using ($10*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'filter' fillstyle pattern 1 lt -1,\
  '' using ($11*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'dual'   fillstyle pattern 2 lt -1,\
  '' using ($12*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'match'  fillstyle pattern 3 lt rgb "#8B0000",\
  '' using ($13*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'kernel' fillstyle pattern 3 lt -1,\
  '' using ($14*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'repair' fillstyle pattern 2 lt rgb "#8B0000"

file = '_4/perf_hsw_shock.dat' 
set output "plots/hsw_shock_repar_4.eps"
plot newhistogram '',\
file using ($10*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'primal' fillstyle pattern 1 lt -1,\
  '' using ($11*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'color'  fillstyle pattern 3 lt rgb "#8B0000",\
  '' using ($12*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'qualit' fillstyle pattern 7 lt -1,\
  '' using ($13*1e2/($10+$11+$12+$13+$14)):xticlabels(1) title 'laplacian' fillstyle pattern 2 lt -1
