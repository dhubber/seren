#Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 "Times New Roman" 14
set output "STATPOLY1-AB-CONSTH.ps"
set size 1.0, 1.0
set origin 0.0,0.0
set pointsize 0.2
unset key
set xlabel 'r'

set origin 0.0,0.0
set xrange[0:1]
set yrange[0:2]
#set title 'Density'
#set label "(a)" at -1.85,1.12
set ylabel "{/Symbol r}"
plot 'STATPOLY1-AB-CONSTH.debug.fin' u 17:11 w d, 'poly1.dat' u 6:7 w l
unset label

set terminal x11
quit

