# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
# ----------------------------------------------------------------------------
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "GEO1-NBODY.ps"
set size 1.2, 0.6
set origin 0.0,0.0
set pointsize 0.6
unset title
#unset key
set multiplot

set size 0.6,0.6
set origin 0.0,0.0
set xrange[0.1:1.0]
set yrange[0.000001:1.0]
set logscale
set label 1 "(a)" at 0.108,0.5
set ylabel 'RMS fractional force error, {/Symbol e}'
set format y "10^{%L}"
set xlabel '{/Symbol q}_{MAC}'
set pointsize 1.2
plot 'GEO1-MONO-NBODY/GEO1-MONO-NBODY.dat' u 1:3 title 'Monopole' w lp ls 12 ,\
 'GEO1-QUAD-NBODY/GEO1-QUAD-NBODY.dat' u 1:3 title 'Quadrupole' w lp ls 7,\
 'GEO1-OCT-NBODY/GEO1-OCT-NBODY.dat' u 1:3 title 'Octupole' w lp ls 10

set size 0.6,0.6
set origin 0.6,0.0
set yrange[0.001:2.0]
set xrange[0.000001:0.1]
set logscale
set label 1 "(b)" at 0.0000015,1.2
set xlabel 'RMS fractional force error, {/Symbol e}_{}'
set format y "10^{%L}"
set format x "10^{%L}"
set ylabel 't_{tree}/t_{direct}'
set pointsize 1.2
plot 'GEO1-MONO-NBODY/GEO1-MONO-NBODY.dat' u 3:2 title 'Monopole' w lp ls 12 ,\
'GEO1-QUAD-NBODY/GEO1-QUAD-NBODY.dat' u 3:2 title 'Quadrupole' w lp ls 7, \
'GEO1-OCT-NBODY/GEO1-OCT-NBODY.dat' u 3:2 title 'Octupole' w lp ls 10

unset multiplot
set terminal x11
quit
