# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "STATPOLY1-AB-GRADH.ps"
set size 0.6, 0.6
set origin 0.0,0.0
set pointsize 0.2
unset key
set xlabel 'r'

set origin 0.0,0.0
set xrange[0:1]
set yrange[0:2]
#set title 'Density'
#set label "(a)" at -1.85,1.12
set label "(b)" at 0.92,1.85
set ylabel "{/Symbol r}"
plot 'STATPOLY1-AB-GRADH.debug.fin' u 17:11 w p ls 7, 'poly1.dat' u 6:7 w l lt 1 lc rgbcolor "red"
unset label

set terminal x11
quit

