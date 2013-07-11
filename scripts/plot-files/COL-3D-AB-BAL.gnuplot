# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "COL-3D-AB-BAL.ps"
set size 1.2, 0.6
set origin 0.0,0.0
set pointsize 0.2
unset key
set multiplot 
set xlabel 'x'

set size 0.6,0.6
set origin 0.0,0.0
set xrange[-0.5:0.5]
set yrange[0.0:20.0]
set label "(a)" at 0.4,18.0
set ylabel "{/Symbol r}"
plot 'COL-3D-AB-BAL.debug.fin' u 1:11 w p 7, 'COL.exact' u 1:2 w l 1
unset label

set size 0.6,0.6
set origin 0.6,0.0
set ylabel "v_x"
set yrange[-5.0:5.0]
set label "(b)" at 0.4,4.0
plot 'COL-3D-AB-BAL.debug.fin' u 1:4 w p 7, 'COL.exact' u 1:3 w l 1
unset label

unset multiplot
set terminal x11
quit

