# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "ISOFREEFALL1-GRADH-BH.ps"
#set size 1.2, 0.6
set size 0.6, 0.6
set origin 0.0,0.0
set pointsize 0.6
unset key
#set multiplot

set size 0.6, 0.6
set origin 0.0,0.0
set xlabel 't/t_{ff}'
set xrange[0:1]
set yrange[0:1.2]
set ylabel 'r/R_0'
#set label "(a)" at 0.9,1.1
set label "10%" at 0.02,0.52
set label "50%" at 0.02,0.85
set label "90%" at 0.02,1.06
plot 'ISOFREEFALL1-GRADH-BH.gasradius' every 1 u ($1/1.1107207):($10) w l lt 3, \
'freefall.dat' u ($1):(0.965489*$2) w l lt 1, \
'ISOFREEFALL1-GRADH-BH.gasradius' every 1 u ($1/1.1107207):($6) w l lt 3, \
'freefall.dat' u ($1):(0.7937005*$2) w l lt 1, \
'ISOFREEFALL1-GRADH-BH.gasradiud' every 1 u ($1/1.1107207):($2) w l lt 3, \
'freefall.dat' u ($1):(0.4641588*$2) w l lt 1, \
'isofreefall.dat' u 1:2 w l lt 6 lc rgbcolor "black"
unset label

#set size 0.6, 0.6
#set origin 0.6, 0.0
#set xlabel 'r/R_0'
#set ylabel '{/Symbol r}/{/Symbol r}_0'
#set xrange [0.01:2.0]
#set yrange[0.2:100]
#set label "(b)" at 1.2,60
##set label "t = 0" at 0.05,1.05
##set label "t = 0.7t_{ff}" at 0.05,4.0
##set label "t = 0.9t_{ff}" at 0.05,26.0
#set logscale
#plot 'ISOFREEFALL1-GRADH-BH.debug.ini' u ($17):($11/0.238732414) w d, 'ISOFREEFALL1-GRADH-BH.debug.00014' u ($17):($11/0.238732414) w d, 'ISOFREEFALL1-GRADH-BH.debug.00018' u ($17):($11/0.238732414) w d

#unset multiplot
set terminal x11
quit

