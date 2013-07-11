# ----------------------------------------------------------------------------
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "EIGEN1-NBODY.ps"
set size 1.2, 0.6
set origin 0.0,0.0
set pointsize 0.6
#unset key
#unset title
set multiplot

#unset label
set size 0.6,0.6
set origin 0.0,0.0
set xrange[0.000001:0.001]
set yrange[0.000001:1.0]
set logscale
set label 1 "(a)" at 0.00000125,0.5
set ylabel 'RMS fractional force error, {/Symbol e}'
set format y "10^{%L}"
set format x "10^{%L}"
set xlabel '{/Symbol a}_{MAC}'
set pointsize 1.2
plot 'EIGEN1-QUAD-NBODY/EIGEN1-QUAD-NBODY.dat' u 1:3 title 'Quadrupole' w lp ls 7,  'EIGEN1-OCT-NBODY/EIGEN1-OCT-NBODY.dat' u 1:3 title 'Octupole' w lp ls 10

#unset label
set autoscale

set size 0.6,0.6
set origin 0.6,0.0
set xrange[0.000001:0.1]
set yrange[0.001:2.0]
set logscale
set label 1 "(b)" at 0.0000015,1.2
set xlabel 'RMS fractional force error, {/Symbol e}_{}'
set format y "10^{%L}"
set format x "10^{%L}"
set ylabel 't_{tree}/t_{direct}'
set pointsize 1.2
plot 'EIGEN1-QUAD-NBODY/EIGEN1-QUAD-NBODY.dat' u 3:2 title 'Quadrupole' w lp ls 7,  'EIGEN1-OCT-NBODY/EIGEN1-OCT-NBODY.dat' u 3:2 title 'Octupole' w lp ls 10
unset multiplot
set terminal x11
quit
