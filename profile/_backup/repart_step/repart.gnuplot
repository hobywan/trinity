reset
#set terminal postscript eps enhanced color 14 size 6.5cm,6.5cm

kernel=3
set terminal postscript eps enhanced color 14 size 8cm,6.5cm
set output "repart_".kernel.".eps"
set title "Time ratio per step [Kernel ".kernel."]"
# =============================================
set size 1,1
set xlabel ""
set ylabel ""
set yrange [0:100]
   # set xtics rotate by -45
#set ytics 20
#set key outside right box
#set key outside bottom box
# =============================================
set style data histogram
unset key
#set style histogram rowstack gap 1
set style fill solid border -1
set style fill solid 1
#set boxwidth 0.25 relative
set boxwidth 0.25
set style line 2 lc rgb 'black' lt 1 lw 1
set grid y
set style histogram rowstacked
set xtics rotate by 45 right
#set style fill pattern border -1
#iset ytics 10 nomirror
#set logscale y
file = '../benchs/repart_'.kernel.'.dat' 
plot newhistogram '',\
file using 2:xticlabels(1) title 'step 1' fillstyle pattern 1  lt -1,\
'' using 3:xticlabels(1) title 'step 2' fillstyle pattern 3  lt -1,\
'' using 4:xticlabels(1) title 'step 3' fillstyle pattern 7  lt -1,\
'' using 5:xticlabels(1) title 'step 4' fillstyle pattern 10 lt -1  

# =============================================
reset
set terminal postscript eps enhanced color 14 size 8cm,5.5cm
set output "repart_knl.eps"
set title "Ratio time per step - Intel KNL"
# =============================================
set size 1,1
set xlabel ""
set ylabel ""
set yrange [0:100]
   # set xtics rotate by -45
#set ytics 20
#set key outside right box
#set key outside bottom box
# =============================================
set style data histogram
unset key
#set style histogram rowstack gap 1
set style fill solid border -1
set style fill solid 1
#set boxwidth 0.25 relative
set boxwidth 0.25
set style line 2 lc rgb 'black' lt 1 lw 1
set grid y
set style histogram rowstacked
#set xtics rotate by 45 right
#set style fill pattern border -1
#iset ytics 10 nomirror
#set logscale y
file = 'repart_knl.dat' 
plot newhistogram '',\
file using 2:xticlabels(1) title 'step 1' fillstyle pattern 1  lt -1,\
'' using 3:xticlabels(1) title 'step 2' fillstyle pattern 3  lt -1,\
'' using 4:xticlabels(1) title 'step 3' fillstyle pattern 7  lt -1,\
'' using 5:xticlabels(1) title 'step 4' fillstyle pattern 10 lt -1  

# =============================================
reset
set terminal postscript eps enhanced color 14 size 8cm,5.5cm
set output "repart_hsw.eps"
set title "Ratio time per step - Intel Haswell"
# =============================================
set size 1,1
set xlabel ""
set ylabel ""
set yrange [0:100]
   # set xtics rotate by -45
#set ytics 20
#set key outside right box
#set key outside bottom box
# =============================================
set style data histogram
unset key
#set style histogram rowstack gap 1
set style fill solid border -1
set style fill solid 1
#set boxwidth 0.25 relative
set boxwidth 0.25
set style line 2 lc rgb 'black' lt 1 lw 1
set grid y
set style histogram rowstacked
#set xtics rotate by 45 right
#set style fill pattern border -1
#iset ytics 10 nomirror
#set logscale y
file = 'repart_hsw.dat' 
plot newhistogram '',\
file using 2:xticlabels(1) title 'step 1' fillstyle pattern 1  lt -1,\
'' using 3:xticlabels(1) title 'step 2' fillstyle pattern 3  lt -1,\
'' using 4:xticlabels(1) title 'step 3' fillstyle pattern 7  lt -1,\
'' using 5:xticlabels(1) title 'step 4' fillstyle pattern 10 lt -1  	       
