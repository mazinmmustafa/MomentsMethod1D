set terminal wxt
set xrange [0:180]
set yrange [-30:0]
set xlabel("theta [deg]")
set ylabel("RCS [dB]")
plot "DataRCS.dat" using 1:3 with lines lw 2 lc 16 