# Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "SHEAR-3D-GRADH-AB.ps"
set size 0.8,0.8
set origin 0.0,0.0
set pointsize 0.2
unset key

#set origin 0.5,0.0
set xrange[0.0:1.0]
set yrange[-1.2:1.2]
#set title 'Specific energy'
set xlabel 'y'
set ylabel "v_x"
#set label "(a)" at 0.015,3.75
plot 'SHEAR-3D-GRADH-AB.debug.fin' u 2:4 w p 7, sin(6.2831853*x) w l 1
unset label

set terminal x11
quit

