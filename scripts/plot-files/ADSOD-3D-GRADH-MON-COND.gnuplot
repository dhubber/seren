# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "ADSOD-3D-GRADH-MON-COND.ps"
set size 1.2, 1.2
set origin 0.0,0.0
set pointsize 0.2
unset key
set multiplot 
set xlabel 'x'

set size 0.6,0.6
set origin 0.0,0.6
set xrange[-2:2]
set yrange[0.2:1.2]
#set title 'Density'
set label "(a)" at -1.85,1.12
set ylabel "{/Symbol r}"
plot 'ADSOD-3D-GRADH-MON-COND.debug.fin' every 5 u 1:11 w p ls 7, 'ADSOD.exact' u 1:2 w l lt 1 lc rgbcolor "red"
unset label

set size 0.6,0.6
set origin 0.6,0.6
set yrange[-0.1:0.9]
#set title 'x-velocity'
set ylabel "v_x"
set label "(b)" at -1.85,0.82
plot 'ADSOD-3D-GRADH-MON-COND.debug.fin' every 5 u 1:4 w p ls 7, 'ADSOD.exact' u 1:3 w l lt 1 lc rgbcolor "red"
unset label

set size 0.6,0.6
set origin 0.0,0.0
set yrange[0.0:1.2]
#set title 'Pressure'
set ylabel "P"
set label "(c)" at -1.85,1.12
plot 'ADSOD-3D-GRADH-MON-COND.debug.fin' every 5 u 1:16 w p ls 7, 'ADSOD.exact' u 1:4 w l lt 1 lc rgbcolor "red"
unset label

set size 0.6,0.6
set origin 0.6,0.0
set yrange[1.7:2.7]
#set title 'Specific internal energy'
set ylabel "u"
set label "(d)" at -1.85,2.62
plot 'ADSOD-3D-GRADH-MON-COND.debug.fin' every 5 u 1:17 w p ls 7, 'ADSOD.exact' u 1:5 w l lt 1 lc rgbcolor "red"
unset label

unset multiplot
set terminal x11
quit

