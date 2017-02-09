
kernel=3

reset
set terminal postscript eps enhanced color 14 size 8cm,6.5cm
set output "imbalance_".kernel.".eps"
#set size ratio 1
set size 1,1
unset logscale x
set key Left reverse left
set title "Ratio of imbalance [Swapping]"
set ylabel "(%)"
set xlabel "cores"
set yrange [0:50]
#set logscale y 5
#set ytics 8

set style data histogram
set style fill pattern border
set style histogram clustered
set grid
# ---------------------------------
plot "../benchs/nehalem_".kernel.".dat"   using (100*(($5-$6)/$5)):xticlabels(1) title "Intel Nehalem"      lt -1,\
     "../benchs/hasw_nor_".kernel.".dat"  using (100*(($5-$6)/$5)):xticlabels(1) title "Intel Haswell"      lt -1,\
     "../benchs/hasw_hyp_".kernel.".dat"  using (100*(($5-$6)/$5)):xticlabels(1) title "Intel Haswell hyper"lt -1,\
     "../benchs/shock_knl_".kernel.".dat" using (100*(($5-$6)/$5)):xticlabels(1) title "Intel KNL MCD"      lt -1,\
     "../benchs/shock_llc_".kernel.".dat" using (100*(($5-$6)/$5)):xticlabels(1) title "Intel KNL DDR4"     lt -1,\
     "../benchs/shock_hyp_".kernel.".dat" using (100*(($5-$6)/$5)):xticlabels(1) title "Intel KNL hyper"    lt -1


# ---------------------------------
