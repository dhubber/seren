#Graph 1 for sink correction forces paper
#set size 1.2, 0.6
set terminal postscript enhanced color dashed lw 1 #"Times New Roman" 14
set output "SPITZER3-HEALPIX-GRADH.ps"
set size 1.0, 1.0
set origin 0.0,0.0
set pointsize 0.2
#unset key
set xlabel 'r'

set origin 0.0,0.0
set xrange[0:0.12]
set yrange[0.1:1.2]
set title 'R (pc)'
#set label "(a)" at -1.85,1.12
set ylabel "time (Myr)"
plot 'SPITZER3-HEALPIX-GRADH.particles' u 2:5 title 'SPH simulation' w l ls 7 , 0.189*(1.0 + 7.0*11.0*x/4.0/0.189/0.9778)**(4.0/7.0) title 'Spitzer solution' ls 1
unset label

set terminal x11
quit

