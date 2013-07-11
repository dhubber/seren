# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "POLYRAD1-AB.ps"
set size 1.0, 1.0
set origin 0.0,0.0
set pointsize 0.2
unset key
set xlabel 'log ({/Symbol r}/g cm^{-3})'

set origin 0.0,0.0
set xrange[-17:0]
set yrange[0:5]
#set title 'Density'
#set label "(a)" at -1.85,1.12
set ylabel "log (T/K)"
plot 'POLYRAD1-AB.rad' u 7:8 w l lt 1 lc rgbcolor "black", 'mmi.denstemp.dat' u 1:2 w l lt 1 lc rgbcolor "red"
unset label

set terminal x11
quit

