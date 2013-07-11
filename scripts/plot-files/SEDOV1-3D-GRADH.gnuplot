# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "SEDOV1-3D-GRADH.ps"
set size 1.0, 1.0
set origin 0.0,0.0
set pointsize 0.2
unset key

set size 0.6,0.6
#set origin 0.5,0.0
set xrange[0.0:0.4]
set yrange[0.0:4.2]
#set title 'Specific energy'
set xlabel 'x'
set ylabel "{/Symbol r}"
set label "(a)" at 0.015,3.9
plot 'SEDOV1-3D-GRADH.debug.fin' every 5 using 19:11 with dots ls 7, 'SEDOV.exact' u 1:3 w l lt 1 lc rgbcolor "red"
unset label

set terminal x11
quit

