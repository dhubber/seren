# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "FREEFALL1-GRADH-BH.ps"
set size 1.2, 0.6
set origin 0.0,0.0
set pointsize 0.6
unset key
set multiplot

set size 0.6, 0.6
set origin 0.0,0.0
set xlabel 't/t_{ff}'
set xrange[0:1]
set yrange[0:1.2]
set ylabel 'r/R_0'
set label "(a)" at 0.9,1.1
set label "10%" at 0.03,0.52
set label "50%" at 0.03,0.85
set label "90%" at 0.03,1.02
plot 'freefall.dat' u ($1):(0.965489*$2) w l lt 1,\
'FREEFALL1-GRADH-BH.freefall' every 1 u 1:91 w p lc rgbcolor "black" pt 12,\
'freefall.dat' u ($1):(0.7937005*$2) w l lt 1,\
'FREEFALL1-GRADH-BH.freefall' every 1 u 1:51 w p lc rgbcolor "black" pt 12,\
'freefall.dat' u ($1):(0.4641588*$2) w l lt 1,\
'FREEFALL1-GRADH-BH.freefall' every 1 u 1:11 w p lc rgbcolor "black" pt 12
unset label

set size 0.6, 0.6
set origin 0.6,0.0
set xlabel '1 - t/t_{ff}'
set xrange[0.001:1] reverse
set yrange[0.01:1.0]
#set autoscale
set logscale
set ylabel 'r/R_0'
set label "(b)" at 0.002,0.7
set label "10%" at 0.31,0.2
set label "50%" at 0.31,0.35
set label "90%" at 0.31,0.79
plot 'freefall.dat' u (1 - $1):(0.965489*$2) w l lt 1,\
'FREEFALL1-GRADH-BH.freefall' every 1 u (1 - $1):($91) w p lc rgbcolor "black" pt 12,\
'freefall.dat' u (1 - $1):(0.7937005*$2) w l lt 1,\
'FREEFALL1-GRADH-BH.freefall' every 1 u (1 - $1):($51) w p lc rgbcolor "black" pt 12,\
'freefall.dat' u (1 - $1):(0.4641588*$2) w l lt 1,\
'FREEFALL1-GRADH-BH.freefall' every 1 u (1 - $1):($11) w p lc rgbcolor "black" pt 12
unset label


#set size 0.6, 0.6
#set origin 0.6,0.0
#set xlabel 't/t_{ff}'
#set xrange[0.95:1]
#set yrange[0:0.3]
#set ylabel 'r/R_0'
#set label "(b)" at 0.99,0.27
##set label "10%" at 0.025,0.6
##set label "50%" at 0.025,0.9
##set label "90%" at 0.025,1.02
#plot 'freefall.dat' u ($1):(0.965489*$2) w l lt 1,\
#'FREEFALL1-GRADH-BH.freefall' u 1:91 w p lc "black" pt 12,\
#'freefall.dat' u ($1):(0.7937005*$2) w l lt 1,\
#'FREEFALL1-GRADH-BH.freefall' u 1:51 w p lc "black" pt 12,\
#'freefall.dat' u ($1):(0.4641588*$2) w l lt 1,\
#'FREEFALL1-GRADH-BH.freefall' u 1:11 w p lc "black" pt 12
#unset label


unset multiplot
set terminal x11
quit

